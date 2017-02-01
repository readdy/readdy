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
 * @file SingleCPUKernelStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 * @todo
 */

#include <algorithm>
#include <readdy/model/Particle.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy {
namespace kernel {
namespace scpu {
struct SCPUStateModel::Impl {
    double currentEnergy = 0;
    std::unique_ptr<model::SCPUParticleData> particleData;
    std::unique_ptr<model::SCPUNeighborList> neighborList;
    std::vector<std::unique_ptr<readdy::model::top::Topology>> topologies;
    readdy::model::KernelContext const *context;
};

SCPUStateModel::SCPUStateModel(readdy::model::KernelContext const *context) : pimpl(
        std::make_unique<SCPUStateModel::Impl>()) {
    pimpl->particleData = std::make_unique<model::SCPUParticleData>();
    pimpl->neighborList = std::make_unique<model::SCPUNeighborList>(context);
    pimpl->context = context;
}

void readdy::kernel::scpu::SCPUStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->particleData->addParticles({p});
}

void readdy::kernel::scpu::SCPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->particleData->addParticles(p);
}

const std::vector<readdy::model::Vec3>
readdy::kernel::scpu::SCPUStateModel::getParticlePositions() const {
    const auto &data = *pimpl->particleData;
    std::vector<readdy::model::Vec3> target{};
    target.reserve(data.size());
    for (const auto &entry : data) {
        if(!entry.is_deactivated()) target.push_back(entry.position());
    }
    return target;
}

model::SCPUParticleData *readdy::kernel::scpu::SCPUStateModel::getParticleData() const {
    return pimpl->particleData.get();
}

void readdy::kernel::scpu::SCPUStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->particleData->removeParticle(p);
}

double readdy::kernel::scpu::SCPUStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

void readdy::kernel::scpu::SCPUStateModel::increaseEnergy(double increase) {
    pimpl->currentEnergy += increase;
}

const model::SCPUNeighborList *SCPUStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

const std::vector<readdy::model::Particle> SCPUStateModel::getParticles() const {
    const auto &data = *pimpl->particleData;
    std::vector<readdy::model::Particle> result;
    result.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.is_deactivated()) {
            result.push_back(data.toParticle(entry));
        }
    }
    return result;
}

void SCPUStateModel::updateNeighborList() {
    pimpl->neighborList->create(*pimpl->particleData);
}

void SCPUStateModel::calculateForces() {
    pimpl->currentEnergy = 0;
    // update forces and energy order 1 potentials
    {
        const readdy::model::Vec3 zero{0, 0, 0};
        for(auto &e : *pimpl->particleData) {
            e.force = zero;
            for (const auto &po1 : pimpl->context->getOrder1Potentials(e.type)) {
                po1->calculateForceAndEnergy(e.force, pimpl->currentEnergy, e.position());
            }
        }
    }

    // update forces and energy order 2 potentials
    {
        const auto &difference = pimpl->context->getShortestDifferenceFun();
        readdy::model::Vec3 forceVec{0, 0, 0};
        for (auto &&it = pimpl->neighborList->begin(); it != pimpl->neighborList->end(); ++it) {
            auto i = it->idx1;
            auto j = it->idx2;
            auto& entry_i = pimpl->particleData->entry_at(i);
            auto& entry_j = pimpl->particleData->entry_at(j);
            const auto &potentials = pimpl->context->getOrder2Potentials(entry_i.type, entry_j.type);
            for (const auto &potential : potentials) {
                potential->calculateForceAndEnergy(forceVec, pimpl->currentEnergy, difference(entry_i.position(), entry_j.position()));
                entry_i.force += forceVec;
                entry_j.force += -1 * forceVec;
            }
        }
    }
    // update forces and energy for topologies
    {
        for(const auto& topology : pimpl->topologies) {
            // calculate bonded potentials
            for(const auto& bondedPot : topology->getBondedPotentials()) {
                // auto energy = bondedPot->createForceAndEnergyAction(pimpl->kernel)->perform();
                // pimpl->currentEnergy += energy;
            }
        }
    }
}

void SCPUStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void SCPUStateModel::removeAllParticles() {
    pimpl->particleData->clear();
}

readdy::model::top::Topology *const SCPUStateModel::addTopology(std::vector<readdy::model::Particle> &particles) {
    for(const auto& p : particles) {
        if(pimpl->context->getParticleTypeInfo(p.getType()).flavor != readdy::model::Particle::FLAVOR_TOPOLOGY) {
            throw std::invalid_argument("at least one of the particles did not have the topology flavor!");
        }
    }
    std::vector<std::size_t> ids = pimpl->particleData->addParticlesGetIds(particles);
    pimpl->topologies.push_back(std::make_unique<readdy::model::top::Topology>(std::move(ids)));
    return pimpl->topologies.back().get();
}


SCPUStateModel &SCPUStateModel::operator=(SCPUStateModel &&rhs) = default;

SCPUStateModel::SCPUStateModel(SCPUStateModel &&rhs) = default;

SCPUStateModel::~SCPUStateModel() = default;
}
}
}




