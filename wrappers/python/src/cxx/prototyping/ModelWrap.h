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

#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace py = pybind11;

namespace readdy {
namespace rpy {
class Model : public readdy::kernel::scpu::SCPUStateModel {

    using super = readdy::kernel::scpu::SCPUStateModel;
public:

    using super::SCPUStateModel;

    virtual void removeParticle(const readdy::model::Particle &p) override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(void, super, "remove_particle", removeParticle, p);
    }

    virtual const std::vector<model::Vec3> getParticlePositions() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(const std::vector<model::Vec3>, super, "get_particle_positions", getParticlePositions,);
    }

    virtual double getEnergy() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(double, super, "get_energy", getEnergy,);
    }

    virtual void increaseEnergy(double increase) override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(void, super, "increase_energy", increaseEnergy, increase);
    }

    virtual kernel::scpu::model::SCPUParticleData *getParticleData() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(kernel::scpu::model::SCPUParticleData*, super, "get_particle_data",
                               getParticleData,);
    }

    virtual const readdy::kernel::scpu::model::SCPUNeighborList *getNeighborList() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(const kernel::scpu::model::SCPUNeighborList*, super, "get_neighbor_list",
                               getNeighborList,);
    }

    virtual const std::vector<model::Particle> getParticles() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_NAME(const std::vector<model::Particle>, super, "get_particles", getParticles,);
    }
};
}
}

#endif //READDY_MAIN_MODELWRAP_H
