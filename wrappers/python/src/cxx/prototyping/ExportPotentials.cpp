//
// Created by mho on 10/08/16.
//
#include <boost/python.hpp>
#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/potentials/PotentialFactory.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder1.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder2.h>
#include "../PyConverters.h"

namespace bpy = boost::python;
namespace mpl = boost::mpl;
namespace rpy = readdy::py;
namespace pot = readdy::model::potentials;
namespace spot = readdy::kernel::singlecpu::potentials;

using _rdy_vec = readdy::model::Vec3;

using _rdy_pot = pot::Potential;
using _rdy_pot1 = pot::PotentialOrder1;
using _rdy_pot2 = pot::PotentialOrder2;

using _rdy_pot_factory = pot::PotentialFactory;

struct PyPotentialO1 : public _rdy_pot1, bpy::wrapper<_rdy_pot1> {
    PyPotentialO1(std::string name) : _rdy_pot1(name) {}

    virtual double calculateEnergy(const _rdy_vec &position) const override {
        return this->get_override("calculate_energy")(position);
    }

    _rdy_vec internalCalculateForce(const _rdy_vec &position) const {
        return this->get_override("calculate_force")(position);
    }

    virtual void calculateForce(_rdy_vec &force, const _rdy_vec &position) const override {
        force += internalCalculateForce(position);
    }

    virtual void calculateForceAndEnergy(_rdy_vec &force, double &energy,
                                         const _rdy_vec &position) const override {
        calculateForce(force, position);
        energy += calculateEnergy(position);
    }

    virtual void configureForType(const unsigned int type) override {
        this->get_override("configure_for_type")(type);
    }

    virtual double getRelevantLengthScale() const noexcept override {
        return this->get_override("get_relevant_length_scale")();
    }

    virtual double getMaximalForce(double kbt) const noexcept override {
        return this->get_override("get_maximal_force")(kbt);
    }

    virtual PyPotentialO1 *replicate() const override {
        return new PyPotentialO1(*this);
    }
};

struct PyPotentialO2 : public _rdy_pot2, bpy::wrapper<_rdy_pot2> {
    PyPotentialO2(std::string name) : PotentialOrder2(name) {}

    virtual double getMaximalForce(double kbt) const noexcept override {
        return this->get_override("get_maximal_force")(kbt);
    }

    virtual double calculateEnergy(const _rdy_vec &x_ij) const override {
        return this->get_override("calculate_energy")(x_ij);
    }

    _rdy_vec internalCalculateForce(const _rdy_vec &x_ij) const {
        return this->get_override("calculate_force")(x_ij);
    }

    virtual void calculateForce(_rdy_vec &force, const _rdy_vec &x_ij) const override {
        force += internalCalculateForce(x_ij);
    }

    virtual void calculateForceAndEnergy(_rdy_vec &force, double &energy, const _rdy_vec &x_ij) const override {
        force += internalCalculateForce(x_ij);
        energy += calculateEnergy(x_ij);
    }

    virtual double getCutoffRadiusSquared() const override {
        const auto cutoff = getCutoffRadius();
        return cutoff * cutoff;
    }

    virtual void configureForTypes(unsigned int type1, unsigned int type2) override {
        this->get_override("configure_for_types")(type1, type2);
    }

    virtual pot::PotentialOrder2 *replicate() const override {
        return new PyPotentialO2(*this);
    }

    virtual double getCutoffRadius() const override {
        return this->get_override("get_cutoff_radius")();
    }
};

template<std::size_t owner = 1, class bp = bpy::default_call_policies>
using internal_ref = bpy::return_internal_reference<owner, bp>;

void exportPotentials() {

    bpy::class_<spot::CubePotential, boost::noncopyable>("CubePotential", bpy::no_init)
            .def("get_name", &spot::CubePotential::getName, bpy::return_value_policy<bpy::copy_const_reference>());
    bpy::class_<spot::HarmonicRepulsion>("HarmonicRepulsion", bpy::no_init)
            .def("get_name", &spot::HarmonicRepulsion::getName, bpy::return_value_policy<bpy::copy_const_reference>());
    bpy::class_<spot::WeakInteractionPiecewiseHarmonic>("WeakInteractionPiecewiseHarmonic", bpy::no_init)
            .def("get_name", &spot::WeakInteractionPiecewiseHarmonic::getName, bpy::return_value_policy<bpy::copy_const_reference>());

    bpy::class_<PyPotentialO1, boost::noncopyable>("PotentialOrder1", bpy::init<std::string>())
            .def("calculate_energy", bpy::pure_virtual(&_rdy_pot1::calculateEnergy))
            .def("calculate_force", bpy::pure_virtual(&PyPotentialO1::internalCalculateForce))
            .def("configure_for_type", bpy::pure_virtual(&_rdy_pot1::configureForType))
            .def("get_relevant_length_scale", bpy::pure_virtual(&_rdy_pot1::getRelevantLengthScale))
            .def("get_maximal_force", bpy::pure_virtual(&_rdy_pot1::getMaximalForce));
    bpy::class_<PyPotentialO2, boost::noncopyable>("PotentialOrder2", bpy::init<std::string>())
            .def("calculate_energy", bpy::pure_virtual(&_rdy_pot2::calculateEnergy))
            .def("calculate_force", bpy::pure_virtual(&PyPotentialO2::internalCalculateForce))
            .def("configure_for_types", bpy::pure_virtual(&_rdy_pot2::configureForTypes))
            .def("get_cutoff_radius", bpy::pure_virtual(&_rdy_pot2::getCutoffRadius))
            .def("get_maximal_force", bpy::pure_virtual(&_rdy_pot2::getMaximalForce));
    
    auto f_create_cube_pot = &_rdy_pot_factory::createPotential<spot::CubePotential>;
    auto f_create_harmonic_pot = &_rdy_pot_factory::createPotential<spot::HarmonicRepulsion>;
    auto f_create_weak_inter_pot = &_rdy_pot_factory::createPotential<spot::WeakInteractionPiecewiseHarmonic>;
    bpy::class_<_rdy_pot_factory>("PotentialFactory", bpy::no_init)
            .def("create_cube_potential", rpy::adapt_unique(f_create_cube_pot))
            .def("create_harmonic_repulsion", rpy::adapt_unique(f_create_harmonic_pot))
            .def("create_weak_interaction", rpy::adapt_unique(f_create_weak_inter_pot));
}