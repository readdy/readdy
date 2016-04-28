/**
 * << detailed description >>
 *
 * @file SingleCPUKernelStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 * @todo
 */

#include <boost/make_unique.hpp>
#include <readdy/model/Particle.h>
#include <vector>
#include <algorithm>
#include "SingleCPUKernelStateModel.h"

namespace k = readdy::kernel::singlecpu;

struct k::SingleCPUKernelStateModel::Impl {
    readdy::model::time_step_type t;
};

void k::SingleCPUKernelStateModel::updateModel(model::time_step_type t, bool forces, bool distances) {
    pimpl->t = t;
    //TODO
}


k::SingleCPUKernelStateModel::SingleCPUKernelStateModel() : pimpl(boost::make_unique<k::SingleCPUKernelStateModel::Impl>()) {
    particleData = std::make_shared<ParticleData>();
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticle(const model::Particle &p) {
    particleData->addParticles({p});
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticles(const std::vector<model::Particle> &p) {
    particleData->addParticles(p);
}

std::vector<readdy::model::Vec3> readdy::kernel::singlecpu::SingleCPUKernelStateModel::getParticlePositions() {
    return {*particleData->positions};
}

readdy::model::time_step_type readdy::kernel::singlecpu::SingleCPUKernelStateModel::getCurrentTimeStep() {
    return pimpl->t;
}


k::SingleCPUKernelStateModel &k::SingleCPUKernelStateModel::operator=(k::SingleCPUKernelStateModel &&rhs) = default;

k::SingleCPUKernelStateModel::SingleCPUKernelStateModel(k::SingleCPUKernelStateModel &&rhs) = default;

k::SingleCPUKernelStateModel::~SingleCPUKernelStateModel() = default;


readdy::kernel::singlecpu::ParticleData::ParticleData() {
    ids = std::make_shared<std::vector<boost::uuids::uuid>>();
    positions = std::make_shared<std::vector<readdy::model::Vec3>>();
    type = boost::make_unique<std::vector<unsigned int>>();
}

void readdy::kernel::singlecpu::ParticleData::addParticles(const std::vector<model::Particle> particles) {
    for(auto&& p : particles) {
        ids->push_back(p.id);
        positions->push_back(p.pos);
        type->push_back(p.type);
    }
}



