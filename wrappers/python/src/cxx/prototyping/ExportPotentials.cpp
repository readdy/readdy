//
// Created by mho on 10/08/16.
//
#include <pybind11/pybind11.h>

#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/potentials/PotentialFactory.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder1.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder2.h>

namespace py = pybind11;
namespace pot = readdy::model::potentials;
namespace spot = readdy::kernel::scpu::potentials;

using rvp = py::return_value_policy;

using rdy_vec = readdy::model::Vec3;

using rdy_pot = pot::Potential;
using rdy_pot1 = pot::PotentialOrder1;
using rdy_pot2 = pot::PotentialOrder2;

using rdy_pot_factory = pot::PotentialFactory;

class PyPotentialO1 : public rdy_pot1 {
private:
public:
    using rdy_pot1::PotentialOrder1;

    virtual double calculateEnergy(const rdy_vec &position) const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(double, rdy_pot1, "calculate_energy", calculateEnergy, position);
    }

    virtual rdy_vec calculateForceInternal(const rdy_vec& pos) const {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(rdy_vec, PyPotentialO1, "calculate_force", calculateForce, pos);
    }

    virtual void calculateForce(rdy_vec &force, const rdy_vec &position) const override {
        force += calculateForceInternal(position);
    }

    virtual void calculateForceAndEnergy(rdy_vec &force, double &energy, const rdy_vec &position) const override {
        calculateForce(force, position);
        energy += calculateEnergy(position);
    }

    virtual void configureForType(const unsigned int type) override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(void, PyPotentialO1, "configure_for_type", configureForType, type);
    }

    virtual double getRelevantLengthScale() const noexcept override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(double, rdy_pot1, "get_relevant_length_scale", getRelevantLengthScale,);
    }

    virtual double getMaximalForce(double kbt) const noexcept override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(double, rdy_pot1, "get_maximal_force", getMaximalForce, kbt);
    }

};

class PyPotentialO2 : public rdy_pot2 {
public:

    using rdy_pot2::PotentialOrder2;

    virtual double getMaximalForce(double kbt) const noexcept override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(double, rdy_pot2, "get_maximal_force", getMaximalForce, kbt);
    }

    virtual double calculateEnergy(const rdy_vec &x_ij) const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(double, rdy_pot2, "calculate_energy", calculateEnergy, x_ij);
    }

    virtual rdy_vec calculateForceInternal(const rdy_vec& pos) const {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(rdy_vec, PyPotentialO2, "calculate_force", calculateForce, pos);
    }

    virtual void calculateForce(rdy_vec &force, const rdy_vec &x_ij) const override {
        force += calculateForceInternal(x_ij);
    }

    virtual void calculateForceAndEnergy(rdy_vec &force, double &energy, const rdy_vec &x_ij) const override {
        calculateForce(force, x_ij);
        energy += calculateEnergy(x_ij);
    }

    virtual double getCutoffRadiusSquared() const override {
        const auto cutoff = getCutoffRadius();
        return cutoff * cutoff;
    }

    virtual void configureForTypes(unsigned int type1, unsigned int type2) override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(void, rdy_pot2, "configure_for_types", configureForType, type1, type2);
    }

    virtual double getCutoffRadius() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(double, rdy_pot2, "get_cutoff_radius", getCutoffRadius);
    }
};

void exportPotentials(py::module &proto) {

    py::class_<rdy_pot>(proto, "Potential");
    py::class_<spot::SCPUCubePotential, rdy_pot>(proto, "CubePotential")
            .def("get_name", &spot::SCPUCubePotential::getName);
    py::class_<spot::HarmonicRepulsion, rdy_pot>(proto, "HarmonicRepulsion")
            .def("get_name", &spot::HarmonicRepulsion::getName);
    py::class_<spot::WeakInteractionPiecewiseHarmonic, rdy_pot>(proto, "WeakInteractionPiecewiseHarmonic")
            .def("get_name", &spot::WeakInteractionPiecewiseHarmonic::getName);

    py::class_<rdy_pot1, PyPotentialO1>(proto, "PotentialOrder1")
            .def(py::init<std::string>())
            .def("calculate_energy", &rdy_pot1::calculateEnergy)
            .def("calculate_force", &rdy_pot1::calculateForce)
            .def("configure_for_type", &rdy_pot1::configureForType)
            .def("get_relevant_length_scale", &rdy_pot1::getRelevantLengthScale)
            .def("get_maximal_force", &rdy_pot1::getMaximalForce);
    py::class_<rdy_pot2, PyPotentialO2>(proto, "PotentialOrder2")
            .def(py::init<std::string>())
            .def("calculate_energy", &rdy_pot2::calculateEnergy)
            .def("calculate_force", &rdy_pot2::calculateForce)
            .def("configure_for_types", &rdy_pot2::configureForTypes)
            .def("get_cutoff_radius", &rdy_pot2::getCutoffRadius)
            .def("get_maximal_force", &rdy_pot2::getMaximalForce);

    auto f_create_cube_pot = &rdy_pot_factory::createPotential<spot::SCPUCubePotential>;
    auto f_create_harmonic_pot = &rdy_pot_factory::createPotential<spot::HarmonicRepulsion>;
    auto f_create_weak_inter_pot = &rdy_pot_factory::createPotential<spot::WeakInteractionPiecewiseHarmonic>;
    py::class_<rdy_pot_factory>(proto, "PotentialFactory")
            .def("create_cube_potential", f_create_cube_pot)
            .def("create_harmonic_repulsion", f_create_harmonic_pot)
            .def("create_weak_interaction", f_create_weak_inter_pot);
}