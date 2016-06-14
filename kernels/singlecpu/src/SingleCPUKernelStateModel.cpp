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
#include <readdy/kernel/singlecpu/SingleCPUKernelStateModel.h>
#include <readdy/common/make_unique.h>
#include <boost/log/trivial.hpp>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>
#include <readdy/model/potentials/PotentialOrder2.h>

namespace kern = readdy::kernel::singlecpu;

struct kern::SingleCPUKernelStateModel::Impl {
    readdy::model::time_step_type t = 0;
    double currentEnergy = 0;
    std::unique_ptr<model::SingleCPUParticleData> particleData;
    std::unique_ptr<model::SingleCPUNeighborList> neighborList;
    readdy::model::KernelContext const* context;
    bool firstCall = true;
};

void kern::SingleCPUKernelStateModel::updateModel(readdy::model::time_step_type t, bool forces, bool distances) {
    if(pimpl->firstCall) {
       /* for(auto typePair : pimpl->context->getAllOrder2RegisteredPotentialTypes()) {
            for(auto&& potential : pimpl->context->getOrder2Potentials(std::get<0>(typePair), std::get<1>(typePair))) {
                if(potential->getOrder() == 2) {
                    // todo discuss this: how and when to configure r²
                }
            }
        }*/
    }
    const auto timeStepChanged = t != pimpl->t;
    pimpl->t = t;
    if(timeStepChanged) {
        fireTimeStepChanged();
        pimpl->currentEnergy = 0;
        pimpl->neighborList->create(*pimpl->particleData);
    }

    if(forces) {
        // zero out forces
        const readdy::model::Vec3 zeroVector = {0,0,0};
        std::fill(pimpl->particleData->begin_forces(), pimpl->particleData->end_forces(), zeroVector);
        // update forces and energy
        for(auto&& it = pimpl->neighborList->begin(); it != pimpl->neighborList->end(); ++it) {
            auto i = it->idx1;
            auto j = it->idx2;
            auto type_i = *(pimpl->particleData->begin_types() + i);
            auto type_j = *(pimpl->particleData->begin_types() + j);
            const auto& pos_i = *(pimpl->particleData->begin_positions() + i);
            const auto& pos_j = *(pimpl->particleData->begin_positions() + j);
            const auto&& potentials = pimpl->context->getOrder2Potentials(type_i, type_j);
            for(const auto& potential : potentials) {
                if(potential->getOrder() == 2) {
                    static_cast<readdy::model::potentials::PotentialOrder2 *>(potential)
                            ->calculateForceAndEnergy(
                                    *(pimpl->particleData->begin_forces()+i)
                                    , pimpl->currentEnergy, pos_i, pos_j
                            );
                } else {
                    // TODO
                }
            }
        }
    }

    pimpl->firstCall = false;
}


kern::SingleCPUKernelStateModel::SingleCPUKernelStateModel(readdy::model::KernelContext const* context) : pimpl(std::make_unique<kern::SingleCPUKernelStateModel::Impl>()) {
    pimpl->particleData = std::make_unique<model::SingleCPUParticleData>(10);
    pimpl->neighborList = std::make_unique<model::NaiveSingleCPUNeighborList>();
    pimpl->context = context;
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->particleData->addParticles({p});
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->particleData->addParticles(p);
}

const std::vector<readdy::model::Vec3> readdy::kernel::singlecpu::SingleCPUKernelStateModel::getParticlePositions() const {
    std::vector<readdy::model::Vec3> target {};
    target.reserve(pimpl->particleData->size());
    std::copy(pimpl->particleData->begin_positions(), pimpl->particleData->end_positions(), std::back_inserter(target));
    return target;
}

const readdy::model::time_step_type readdy::kernel::singlecpu::SingleCPUKernelStateModel::getCurrentTimeStep() const {
    return pimpl->t;
}

kern::model::SingleCPUParticleData* readdy::kernel::singlecpu::SingleCPUKernelStateModel::getParticleData() const {
    return pimpl->particleData.get();
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->particleData->removeParticle(p);
}

double readdy::kernel::singlecpu::SingleCPUKernelStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::increaseEnergy(double increase) {
    pimpl->currentEnergy += increase;
}


kern::SingleCPUKernelStateModel &kern::SingleCPUKernelStateModel::operator=(kern::SingleCPUKernelStateModel &&rhs) = default;

kern::SingleCPUKernelStateModel::SingleCPUKernelStateModel(kern::SingleCPUKernelStateModel &&rhs) = default;

kern::SingleCPUKernelStateModel::~SingleCPUKernelStateModel() = default;


