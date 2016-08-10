//
// Created by mho on 10/08/16.
//
#include <boost/python.hpp>
#include <readdy/model/potentials/PotentialOrder1.h>
#include "../PyConverters.h"

namespace bpy = boost::python;
namespace mpl = boost::mpl;
namespace rp = readdy::py;

using _rdy_pot = readdy::model::potentials::Potential;
using _rdy_pot1 = readdy::model::potentials::PotentialOrder1;

struct PyPotentialO1 : public _rdy_pot1, bpy::wrapper<_rdy_pot1> {
    PyPotentialO1(std::string name) : _rdy_pot1(name) {}

    virtual double calculateEnergy(const readdy::model::Vec3 &position) const override {
        return this->get_override("calculate_energy")(position);
    }

    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) const override {
        this->get_override("calculate_force")(force, position);
    }

    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                         const readdy::model::Vec3 &position) const override {
        calculateForce(force, position);
        energy = calculateEnergy(position);
    }

    virtual void configureForType(const unsigned int &type) override {
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

void exportPotentials() {
    bpy::class_<PyPotentialO1, boost::noncopyable>("PotentialOrder1", bpy::init<std::string>())
            .def("calculate_energy", bpy::pure_virtual(&_rdy_pot1::calculateEnergy))
            .def("calculate_force", bpy::pure_virtual(&_rdy_pot1::calculateForce))
            .def("configure_for_type", bpy::pure_virtual(&_rdy_pot1::configureForType))
            .def("get_relevant_length_scale", bpy::pure_virtual(&_rdy_pot1::getRelevantLengthScale))
            .def("get_maximal_force", bpy::pure_virtual(&_rdy_pot1::getMaximalForce));
}