/**
 * << detailed description >>
 *
 * @file CPUStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <future>
#include <readdy/kernel/cpu/CPUStateModel.h>
#include <readdy/common/thread/scoped_thread.h>

namespace readdy {
namespace kernel {
namespace cpu {

namespace thd = readdy::util::thread;

using entries_it = CPUStateModel::data_t::entries_t::iterator;
using neighbors_it = decltype(std::declval<readdy::kernel::cpu::model::NeighborList>().cbegin());
using pot1Map = decltype(std::declval<readdy::model::KernelContext>().getAllOrder1Potentials());
using pot2Map = decltype(std::declval<readdy::model::KernelContext>().getAllOrder2Potentials());
using dist_fun = readdy::model::KernelContext::shortest_dist_fun;

void calculateForcesThread(entries_it begin, entries_it end, neighbors_it neighbors_it,
                           std::promise<double> energyPromise, const CPUStateModel::data_t& data,
                           pot1Map pot1, pot2Map pot2, dist_fun d) {
    double energyUpdate = 0.0;
    for (auto it = begin; it != end; ++it) {
        if (!it->is_deactivated()) {

            //
            // 1st order potentials
            //

            readdy::model::Vec3 force{0, 0, 0};
            const auto& myPos = it->position();

            auto find_it = pot1.find(it->type);
            if(find_it != pot1.end()) {
                for(const auto& potential : find_it->second) {
                    potential->calculateForceAndEnergy(force, energyUpdate, myPos);
                }
            }

            //
            // 2nd order potentials
            //

            for (const auto neighbor : *neighbors_it) {
                auto &neighborEntry = data.entry_at(neighbor);
                auto potit = pot2.find({it->type, neighborEntry.type});
                if (potit != pot2.end()) {
                    auto x_ij = d(myPos, neighborEntry.position());
                    auto distSquared = x_ij * x_ij;
                    for (const auto &potential : potit->second) {
                        if (distSquared < potential->getCutoffRadiusSquared()) {
                            readdy::model::Vec3 updateVec{0, 0, 0};
                            potential->calculateForceAndEnergy(updateVec, energyUpdate, x_ij);
                            force += updateVec;
                        }
                    }
                }
            }

            it->force = force;
        }
        ++neighbors_it;
    }

    energyPromise.set_value(energyUpdate);
}

struct CPUStateModel::Impl {
    readdy::model::KernelContext *context;
    std::unique_ptr<readdy::kernel::cpu::model::NeighborList> neighborList;
    double currentEnergy = 0;

    template<bool fixpos = true>
    const model::ParticleData &cdata() const {
        if (fixpos) particleData->setFixPosFun(context->getFixPositionFun());
        return *particleData;
    }

    template<bool fixpos = true>
    model::ParticleData &data() {
        if (fixpos) particleData->setFixPosFun(context->getFixPositionFun());
        return *particleData;
    }

    Impl(readdy::model::KernelContext *context) {
        particleData = std::make_unique<CPUStateModel::data_t>(context);
        this->context = context;
    }

private:
    std::unique_ptr<readdy::kernel::cpu::model::ParticleData> particleData;
};

void CPUStateModel::calculateForces() {
    pimpl->currentEnergy = 0;
    const auto &particleData = pimpl->cdata<true>();
    const auto potOrder1 = pimpl->context->getAllOrder1Potentials();
    const auto potOrder2 = pimpl->context->getAllOrder2Potentials();
    auto d = pimpl->context->getShortestDifferenceFun();
    {
    std::vector<std::future<double>> energyFutures;
        energyFutures.reserve(config->nThreads());
        {
            std::vector<thd::scoped_thread> threads;
            threads.reserve(config->nThreads());
            const std::size_t grainSize = (pimpl->cdata().size()) / config->nThreads();
            auto it_data_end = pimpl->data<false>().end();
            auto it_data = pimpl->data<false>().begin();
            auto it_nl = pimpl->neighborList->begin();
            for (auto i = 0; i < config->nThreads() - 1; ++i) {
                std::promise<double> energyPromise;
                energyFutures.push_back(energyPromise.get_future());
                threads.push_back(thd::scoped_thread(
                        std::thread(calculateForcesThread, it_data, it_data+grainSize, it_nl, std::move(energyPromise),
                                    std::cref(particleData), potOrder1, potOrder2, d)
                ));
                it_nl += grainSize;
                it_data += grainSize;
            }
            {
                std::promise<double> lastPromise;
                energyFutures.push_back(lastPromise.get_future());
                threads.push_back(thd::scoped_thread(
                        std::thread(calculateForcesThread, it_data, it_data_end, it_nl, std::move(lastPromise),
                                    std::cref(particleData), potOrder1, potOrder2, d)));
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

const std::vector<readdy::model::Vec3> CPUStateModel::getParticlePositions() const {
    const auto &data = pimpl->cdata();
    std::vector<readdy::model::Vec3> target{};
    target.reserve(data.size());
    for (const auto &entry : data) {
        if(!entry.is_deactivated()) target.push_back(entry.position());
    }
    return target;
}

const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
    const auto &data = pimpl->cdata();
    std::vector<readdy::model::Particle> result;
    result.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.is_deactivated()) {
            result.push_back(data.toParticle(entry));
        }
    }
    return result;
}

void CPUStateModel::updateNeighborList() {
    pimpl->neighborList->create();
}

void CPUStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->data<false>().addParticle(p);
}

void CPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->data<false>().addParticles(p);
}

void CPUStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->data<false>().removeParticle(p);
}

double CPUStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

CPUStateModel::CPUStateModel(readdy::model::KernelContext *const context, util::Config const *const config)
        : pimpl(std::make_unique<Impl>(context)), config(config) {
    pimpl->neighborList = std::make_unique<model::NeighborList>(context, *getParticleData(), config);
}

CPUStateModel::data_t *const CPUStateModel::getParticleData() const {
    return &pimpl->data<false>();
}

model::NeighborList *const CPUStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

void CPUStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void CPUStateModel::removeAllParticles() {
    pimpl->data<false>().clear();
}

CPUStateModel::~CPUStateModel() = default;


}
}
}