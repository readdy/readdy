/**
 * << detailed description >>
 *
 * @file KernelContext.cpp
 * @brief Implementation file of the KernelContext.
 * @author clonker
 * @date 18.04.16
 * @todo make proper reference to KernelContext.h, is kBT really indepdendent of t?
 */

#include <readdy/model/KernelContext.h>
#include <readdy/common/make_unique.h>
#include <unordered_map>
#include <boost/log/trivial.hpp>

using namespace readdy::model;

struct readdy::model::KernelContext::Impl {
    uint typeCounter;
    std::unordered_map<std::string, uint> typeMapping;
    double kBT = 0;
    std::array<double, 3> box_size{};
    std::array<bool, 3> periodic_boundary{};
    std::unordered_map<uint, double> diffusionConstants{};
    double timeStep;
};


double KernelContext::getKBT() const {
    return (*pimpl).kBT;
}

void KernelContext::setKBT(double kBT) {
    (*pimpl).kBT = kBT;
}

void KernelContext::setBoxSize(double dx, double dy, double dz) {
    (*pimpl).box_size = {dx, dy, dz};
}

void KernelContext::setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z) {
    (*pimpl).periodic_boundary = {pb_x, pb_y, pb_z};
}

KernelContext::KernelContext() : pimpl(std::make_unique<KernelContext::Impl>()) { }

std::array<double, 3>& KernelContext::getBoxSize() const {
    return pimpl->box_size;
}

std::array<bool, 3>& KernelContext::getPeriodicBoundary() const {
    return pimpl->periodic_boundary;
}

KernelContext::KernelContext(const KernelContext &rhs) : pimpl(std::make_unique<KernelContext::Impl>(*rhs.pimpl)) { }

KernelContext &KernelContext::operator=(const KernelContext &rhs) {
    *pimpl = *rhs.pimpl;
    return *this;
}

double KernelContext::getDiffusionConstant(const std::string& particleType) const {
    return pimpl->diffusionConstants[pimpl->typeMapping[particleType]];
}

void KernelContext::setDiffusionConstant(const std::string& particleType, double D) {
    bool hasType = false;
    if(pimpl->typeMapping.find(particleType) != pimpl->typeMapping.end()) {
        BOOST_LOG_TRIVIAL(warning) << "diffusion constant for particle type " << particleType << " was already set to " << (getDiffusionConstant(particleType)) << " and is now overwritten.";
        hasType = true;
    }
    uint t_id;
    if(hasType) {
        t_id = pimpl->typeMapping[particleType];
    } else {
        t_id = ++(pimpl->typeCounter);
        pimpl->typeMapping.emplace(particleType, t_id);
    }
    pimpl->diffusionConstants[t_id] = D;
}

double KernelContext::getTimeStep() const {
    return pimpl->timeStep;
}

void KernelContext::setTimeStep(double dt) {
    pimpl->timeStep = dt;
}

unsigned int KernelContext::getParticleTypeID(const std::string& name) const {
    return pimpl->typeMapping[name];
}

double KernelContext::getDiffusionConstant(uint particleType) const {
    return pimpl->diffusionConstants[particleType];
}


KernelContext &KernelContext::operator=(KernelContext &&rhs) = default;

KernelContext::KernelContext(KernelContext &&rhs) = default;

KernelContext::~KernelContext() = default;



