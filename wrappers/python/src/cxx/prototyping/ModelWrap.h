/**
 * << detailed description >>
 *
 * @file ModelWrap.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 08.08.16
 */

#ifndef READDY_MAIN_MODELWRAP_H
#define READDY_MAIN_MODELWRAP_H

#include <boost/python.hpp>
#include <readdy/kernel/singlecpu/SingleCPUKernelStateModel.h>

namespace bpy = boost::python;

namespace readdy {
namespace py {
struct Model
        : public readdy::kernel::singlecpu::SingleCPUKernelStateModel,
          bpy::wrapper<readdy::kernel::singlecpu::SingleCPUKernelStateModel> {

    using super = readdy::kernel::singlecpu::SingleCPUKernelStateModel;

    Model(const model::KernelContext *context) : SingleCPUKernelStateModel(context) {}

    virtual void removeParticle(const readdy::model::Particle &p) override {
        if (auto f = this->get_override("remove_particle")) f(p);
        else super::removeParticle(p);
    }

    virtual const std::vector<model::Vec3> getParticlePositions() const override {
        if (auto f = this->get_override("get_particle_positions")) return f();
        return super::getParticlePositions();
    }

    virtual double getEnergy() const override {
        if (auto f = this->get_override("get_energy")) return f();
        return super::getEnergy();
    }

    virtual void increaseEnergy(double increase) override {
        if (auto f = this->get_override("increase_energy")) f(increase);
        else super::increaseEnergy(increase);
    }

    virtual kernel::singlecpu::model::SingleCPUParticleData *getParticleData() const override {
        if (auto f = this->get_override("get_particle_data")) return f();
        return super::getParticleData();
    }

    virtual const readdy::kernel::singlecpu::model::SingleCPUNeighborList *getNeighborList() const override {
        if (auto f = this->get_override("get_neighbor_list")) return f();
        return super::getNeighborList();
    }

    virtual const std::vector<model::Particle> getParticles() const override {
        if (auto f = this->get_override("get_particles")) return f();
        return super::getParticles();
    }

    virtual void default_removeParticle(const readdy::model::Particle &p) {
        super::removeParticle(p);
    }

    virtual const std::vector<model::Vec3> default_getParticlePositions() const {
        return super::getParticlePositions();
    }

    virtual double default_getEnergy() const {
        return super::getEnergy();
    }

    virtual void default_increaseEnergy(double increase) {
        super::increaseEnergy(increase);
    }

    virtual kernel::singlecpu::model::SingleCPUParticleData *default_getParticleData() const {
        return super::getParticleData();
    }

    virtual const kernel::singlecpu::model::SingleCPUNeighborList *default_getNeighborList() const {
        return super::getNeighborList();
    }

    virtual const std::vector<model::Particle> default_getParticles() const {
        return super::getParticles();
    }
};
}
}

#endif //READDY_MAIN_MODELWRAP_H
