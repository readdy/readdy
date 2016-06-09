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

namespace kern = readdy::kernel::singlecpu;

struct kern::SingleCPUKernelStateModel::Impl {
    readdy::model::time_step_type t = 0;
    double currentEnergy = 0;
    std::unique_ptr<model::SingleCPUParticleData> particleData;
};

void kern::SingleCPUKernelStateModel::updateModel(readdy::model::time_step_type t, bool forces, bool distances) {
    const auto timeStepChanged = t != pimpl->t;
    pimpl->t = t;
    //TODO update
    if(timeStepChanged) {
        fireTimeStepChanged();
        pimpl->currentEnergy = 0;
    }
}


kern::SingleCPUKernelStateModel::SingleCPUKernelStateModel() : pimpl(std::make_unique<kern::SingleCPUKernelStateModel::Impl>()) {
    pimpl->particleData = std::make_unique<model::SingleCPUParticleData>(10);
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


