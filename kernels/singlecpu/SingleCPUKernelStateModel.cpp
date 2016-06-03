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
#include <boost/log/trivial.hpp>

namespace k = readdy::kernel::singlecpu;

/*
namespace {
    struct particle_compare : public std::unary_function<readdy::model::Particle, bool> {

        explicit particle_compare(const readdy::model::Particle &p) : p(p) { }

        bool operator()(const readdy::model::Particle &rhs) {
            return rhs.getId() == p.getId();
        }

        readdy::model::Particle p;
    };
}
 */

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

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::removeParticle(const model::Particle &p) {
    auto &&beginIt = pimpl->particleData->ids->begin();
    auto &&endIt = pimpl->particleData->ids->end();
    auto &&it = std::find(beginIt, endIt, p.getId());
    if(it != endIt) {
        const auto idx = it - beginIt;
        pimpl->particleData->deactivatedParticles->push_back(idx);
    } else {
        BOOST_LOG_TRIVIAL(warning) << "Could not find and thus remove particle";
    }
}


k::SingleCPUKernelStateModel &k::SingleCPUKernelStateModel::operator=(k::SingleCPUKernelStateModel &&rhs) = default;

k::SingleCPUKernelStateModel::SingleCPUKernelStateModel(k::SingleCPUKernelStateModel &&rhs) = default;

k::SingleCPUKernelStateModel::~SingleCPUKernelStateModel() = default;


readdy::kernel::singlecpu::ParticleData::ParticleData() {
    ids = std::make_unique<std::vector<boost::uuids::uuid>>();
    positions = std::make_unique<std::vector<readdy::model::Vec3>>();
    type = std::make_unique<std::vector<unsigned int>>();
    forces = std::make_unique<std::vector<readdy::model::Vec3>>();
    deactivatedParticles = std::make_unique<std::vector<unsigned int>>();
}

void readdy::kernel::singlecpu::ParticleData::addParticles(const std::vector<model::Particle> particles) {
    for(auto&& p : particles) {
        ids->push_back(p.getId());
        positions->push_back(p.getPos());
        type->push_back(p.getType());
        forces->push_back({0, 0, 0});
    }
}



