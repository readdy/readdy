/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file StateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */

#include <future>
#include <readdy/kernel/cpu_dense/CPUDStateModel.h>
#include <readdy/common/thread/scoped_thread.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {

namespace thd = readdy::util::thread;

using entries_it = CPUDStateModel::data_t::entries_t::iterator;
using neighbors_it = decltype(std::declval<readdy::kernel::cpu_dense::model::CPUDNeighborList>().cbegin());
using pot1Map = decltype(std::declval<readdy::model::KernelContext>().getAllOrder1Potentials());
using pot2Map = decltype(std::declval<readdy::model::KernelContext>().getAllOrder2Potentials());

void calculateForcesThread(entries_it begin, entries_it end, neighbors_it neighbors_it,
                           std::promise<double> energyPromise, const CPUDStateModel::data_t &data,
                           pot1Map pot1, pot2Map pot2, bool secondOrder) {
    double energyUpdate = 0.0;
    for (auto it = begin; it != end; ++it, ++neighbors_it) {
        //
        // 1st order potentials
        //

        readdy::model::Vec3 force{0, 0, 0};
        const auto &myPos = it->pos;

        auto find_it = pot1.find(it->type);
        if (find_it != pot1.end()) {
            for (const auto &potential : find_it->second) {
                potential->calculateForceAndEnergy(force, energyUpdate, myPos);
            }
        }

        //
        // 2nd order potentials
        //
        if(secondOrder) {
            for (const auto &neighbor : *neighbors_it) {
                auto &neighborEntry = data.entry_at(neighbor.idx);
                auto potit = pot2.find({it->type, neighborEntry.type});
                if (potit != pot2.end()) {
                    for (const auto &potential : potit->second) {
                        if (neighbor.d2 < potential->getCutoffRadiusSquared()) {
                            readdy::model::Vec3 updateVec{0, 0, 0};
                            potential->calculateForceAndEnergy(updateVec, energyUpdate, neighborEntry.pos - it->pos);
                            force += updateVec;
                        }
                    }
                }
            }
        }

        it->force = force;
    }


    energyPromise.set_value(energyUpdate);
}

struct CPUDStateModel::Impl {
    readdy::model::KernelContext *context;
    std::unique_ptr<readdy::kernel::cpu_dense::model::CPUDNeighborList> neighborList;
    double currentEnergy = 0;

    const model::CPUDParticleData &cdata() const {
        return *particleData;
    }

    model::CPUDParticleData &data() {
        return *particleData;
    }

    Impl(readdy::model::KernelContext *context) {
        particleData = std::make_unique<CPUDStateModel::data_t>(context);
        this->context = context;
    }

private:
    std::unique_ptr<readdy::kernel::cpu_dense::model::CPUDParticleData> particleData;
};

void CPUDStateModel::calculateForces() {
    pimpl->currentEnergy = 0;
    const auto &particleData = pimpl->cdata();
    const auto potOrder1 = pimpl->context->getAllOrder1Potentials();
    const auto potOrder2 = pimpl->context->getAllOrder2Potentials();
    {
        std::vector<std::future<double>> energyFutures;
        energyFutures.reserve(config->nThreads());
        {
            std::vector<thd::scoped_thread> threads;
            threads.reserve(config->nThreads());
            const std::size_t grainSize = (pimpl->cdata().size()) / config->nThreads();
            auto it_data_end = pimpl->data().end();
            auto it_data = pimpl->data().begin();
            auto it_nl = pimpl->neighborList->begin();
            const auto data_size = pimpl->data().size();
            const auto nl_size = std::distance(pimpl->neighborList->begin(), pimpl->neighborList->end());
            const auto needSecondOrderPotentials = pimpl->neighborList->getMaxCutoff() > 0;
            if(needSecondOrderPotentials && data_size != nl_size) {
                log::critical("size data = {}, size nl = {}", data_size, nl_size);
            }
            for (auto i = 0; i < config->nThreads() - 1; ++i) {
                std::promise<double> energyPromise;
                energyFutures.push_back(energyPromise.get_future());
                threads.push_back(thd::scoped_thread(
                        std::thread(calculateForcesThread, it_data, it_data + grainSize, it_nl,
                                    std::move(energyPromise),
                                    std::cref(particleData), potOrder1, potOrder2, needSecondOrderPotentials)
                ));
                it_nl += grainSize;
                it_data += grainSize;
            }
            {
                std::promise<double> lastPromise;
                energyFutures.push_back(lastPromise.get_future());
                threads.push_back(thd::scoped_thread(
                        std::thread(calculateForcesThread, it_data, it_data_end, it_nl, std::move(lastPromise),
                                    std::cref(particleData), potOrder1, potOrder2, needSecondOrderPotentials)));
            }

        }
        for (auto &f : energyFutures) {
            pimpl->currentEnergy += f.get();
        }

    }
    /**
     * We need to take 0.5*energy as we iterate over every particle-neighbor-pair twice.
     */
    pimpl->currentEnergy /= 2.0;
}

const std::vector<readdy::model::Vec3> CPUDStateModel::getParticlePositions() const {
    const auto &data = pimpl->cdata();
    std::vector<readdy::model::Vec3> target{};
    target.reserve(data.size());
    for (const auto &entry : data) {
        target.push_back(entry.pos);
    }
    return target;
}

const std::vector<readdy::model::Particle> CPUDStateModel::getParticles() const {
    const auto &data = pimpl->cdata();
    std::vector<readdy::model::Particle> result;
    result.reserve(data.size());
    for (const auto &entry : data) {
        result.push_back(data.toParticle(entry));
    }
    return result;
}

void CPUDStateModel::updateNeighborList() {
    pimpl->neighborList->create();
}

void CPUDStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->data().addParticle(p);
}

void CPUDStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->data().addParticles(p);
}

void CPUDStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->data().removeParticle(p);
}

double CPUDStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

CPUDStateModel::CPUDStateModel(readdy::model::KernelContext *const context, readdy::util::thread::Config const *const config)
        : pimpl(std::make_unique<Impl>(context)), config(config) {
    pimpl->neighborList = std::make_unique<model::CPUDNeighborList>(context, *getParticleData(), config);
}

CPUDStateModel::data_t *const CPUDStateModel::getParticleData() const {
    return &pimpl->data();
}

model::CPUDNeighborList *const CPUDStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

void CPUDStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void CPUDStateModel::removeAllParticles() {
    pimpl->data().clear();
}

readdy::model::top::Topology *const CPUDStateModel::addTopology(const std::vector<readdy::model::TopologyParticle> &particles) {
    return nullptr;
}

CPUDStateModel::~CPUDStateModel() = default;


}
}
}