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
#include "SingleCPUKernelStateModel.h"

namespace k = readdy::kernel::singlecpu;

struct k::SingleCPUKernelStateModel::Impl {
    std::vector<readdy::model::Particle> particles;
};

void k::SingleCPUKernelStateModel::updateModel(bool forces, bool distances) {
    // todo
}


k::SingleCPUKernelStateModel::SingleCPUKernelStateModel() : pimpl(boost::make_unique<k::SingleCPUKernelStateModel::Impl>()) { }

k::SingleCPUKernelStateModel::SingleCPUKernelStateModel(const k::SingleCPUKernelStateModel &rhs) : pimpl(boost::make_unique<k::SingleCPUKernelStateModel::Impl>(*rhs.pimpl)) { }

k::SingleCPUKernelStateModel &k::SingleCPUKernelStateModel::operator=(const k::SingleCPUKernelStateModel &rhs) {
    *pimpl = *rhs.pimpl;
    return *this;
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticle(const model::Particle &p) {
    (*pimpl).particles.push_back(p);
}

void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticles(const std::vector<model::Particle> &p) {
    pimpl->particles.reserve(pimpl->particles.size() + p.size());
    pimpl->particles.insert(pimpl->particles.end(), p.begin(), p.end());
}

std::vector<readdy::model::Particle> readdy::kernel::singlecpu::SingleCPUKernelStateModel::getParticles() const {
    return pimpl->particles;
}


k::SingleCPUKernelStateModel &k::SingleCPUKernelStateModel::operator=(k::SingleCPUKernelStateModel &&rhs) = default;

k::SingleCPUKernelStateModel::SingleCPUKernelStateModel(k::SingleCPUKernelStateModel &&rhs) = default;

k::SingleCPUKernelStateModel::~SingleCPUKernelStateModel() = default;

