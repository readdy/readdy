/**
 * << detailed description >>
 *
 * @file SingleCPUKernelStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 * @todo
 */

#include <readdy/model/Particle.h>
#include <vector>
#include <algorithm>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy {
namespace kernel {
namespace scpu {
struct SCPUStateModel::Impl {
    double currentEnergy = 0;
    std::unique_ptr<model::SCPUParticleData> particleData;
    std::unique_ptr<model::SCPUNeighborList> neighborList;
    readdy::model::KernelContext const *context;
};

SCPUStateModel::SCPUStateModel(readdy::model::KernelContext const *context) : pimpl(
        std::make_unique<SCPUStateModel::Impl>()) {
    pimpl->particleData = std::make_unique<model::SCPUParticleData>(10, true);
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
    const auto size = pimpl->particleData->size();
    std::vector<readdy::model::Vec3> target{};
    target.reserve(size);
    std::copy(pimpl->particleData->begin_positions(), pimpl->particleData->begin_positions() + size,
              std::back_inserter(target));
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
    std::vector<readdy::model::Particle> result;
    for (auto i = 0; i < pimpl->particleData->size(); ++i) {
        result.push_back((*pimpl->particleData)[i]);
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
        auto &&it_forces = pimpl->particleData->begin_forces();
        auto &&it_types = pimpl->particleData->begin_types();
        auto &&it_pos = pimpl->particleData->begin_positions();
        while (it_forces != pimpl->particleData->end_forces()) {
            *it_forces = zero;
            for (const auto &po1 : pimpl->context->getOrder1Potentials(*it_types)) {
                po1->calculateForceAndEnergy(*it_forces, pimpl->currentEnergy, *it_pos);
            }
            ++it_forces;
            ++it_types;
            ++it_pos;
        }
    }

    // update forces and energy order 2 potentials
    {
        const auto &difference = pimpl->context->getShortestDifferenceFun();
        readdy::model::Vec3 forceVec{0, 0, 0};
        for (auto &&it = pimpl->neighborList->begin(); it != pimpl->neighborList->end(); ++it) {
            auto i = it->idx1;
            auto j = it->idx2;
            auto type_i = *(pimpl->particleData->begin_types() + i);
            auto type_j = *(pimpl->particleData->begin_types() + j);
            const auto &pos_i = *(pimpl->particleData->begin_positions() + i);
            const auto &pos_j = *(pimpl->particleData->begin_positions() + j);
            const auto &potentials = pimpl->context->getOrder2Potentials(type_i, type_j);
            for (const auto &potential : potentials) {
                potential->calculateForceAndEnergy(forceVec, pimpl->currentEnergy, difference(pos_i, pos_j));
                *(pimpl->particleData->begin_forces() + i) += forceVec;
                *(pimpl->particleData->begin_forces() + j) += -1 * forceVec;
            }
        }
    }
}

void SCPUStateModel::setNeighborList(std::unique_ptr<model::SCPUNeighborList> ptr) {
    pimpl->neighborList.swap(ptr);
}

void SCPUStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void SCPUStateModel::removeAllParticles() {
    pimpl->particleData->clear();
}


SCPUStateModel &SCPUStateModel::operator=(SCPUStateModel &&rhs) = default;

SCPUStateModel::SCPUStateModel(SCPUStateModel &&rhs) = default;

SCPUStateModel::~SCPUStateModel() = default;
}
}
}




