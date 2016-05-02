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
#include "SingleCPUKernelStateModel.h"
#include <readdy/common/make_unique.h>

namespace k = readdy::kernel::singlecpu;

struct k::SingleCPUKernelStateModel::Impl {
    readdy::model::time_step_type t = 0;
    std::unique_ptr<ParticleData> particleData;
};

void k::SingleCPUKernelStateModel::updateModel(model::time_step_type t, bool forces, bool distances) {
    const auto timeStepChanged = t != pimpl->t;
    pimpl->t = t;
    //TODO update
    if(timeStepChanged) {
        fireTimeStepChanged();
    }
}


k::SingleCPUKernelStateModel::SingleCPUKernelStateModel() : pimpl(std::make_unique<k::SingleCPUKernelStateModel::Impl>()) {
    pimpl->particleData = std::make_unique<ParticleData>();
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticle(const model::Particle &p) {
    pimpl->particleData->addParticles({p});
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticles(const std::vector<model::Particle> &p) {
    pimpl->particleData->addParticles(p);
}

const std::vector<readdy::model::Vec3> readdy::kernel::singlecpu::SingleCPUKernelStateModel::getParticlePositions() const {
    return {*pimpl->particleData->positions};
}

const readdy::model::time_step_type readdy::kernel::singlecpu::SingleCPUKernelStateModel::getCurrentTimeStep() const {
    return pimpl->t;
}

k::ParticleData* readdy::kernel::singlecpu::SingleCPUKernelStateModel::getParticleData() const {
    return pimpl->particleData.get();
}


k::SingleCPUKernelStateModel &k::SingleCPUKernelStateModel::operator=(k::SingleCPUKernelStateModel &&rhs) = default;

k::SingleCPUKernelStateModel::SingleCPUKernelStateModel(k::SingleCPUKernelStateModel &&rhs) = default;

k::SingleCPUKernelStateModel::~SingleCPUKernelStateModel() = default;


readdy::kernel::singlecpu::ParticleData::ParticleData() {
    ids = std::make_shared<std::vector<boost::uuids::uuid>>();
    positions = std::make_shared<std::vector<readdy::model::Vec3>>();
    type = std::make_unique<std::vector<unsigned int>>();
}

void readdy::kernel::singlecpu::ParticleData::addParticles(const std::vector<model::Particle> particles) {
    for(auto&& p : particles) {
        ids->push_back(p.getId());
        positions->push_back(p.getPos());
        type->push_back(p.getType());
    }
}



