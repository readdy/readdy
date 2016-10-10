/**
 * << detailed description >>
 *
 * @file PotentialWrapper.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.06.16
 */

#include "PyPotential.h"

namespace readdy {

namespace py {
PotentialOrder2Wrapper::PotentialOrder2Wrapper(const std::string &name, pybind11::object o1, pybind11::object o2)
        : PotentialOrder2(name),
          calcEnergyFun(new pybind11::object(o1), [](pybind11::object *o) {
              pybind11::gil_scoped_acquire lock;
              delete o;
          }),
          calcForceFun(new pybind11::object(o2), [](pybind11::object *o) {
              pybind11::gil_scoped_acquire lock;
              delete o;
          }) {};

double PotentialOrder2Wrapper::calculateEnergy(const model::Vec3 &x_ij) const {
    pybind11::gil_scoped_acquire lock;
    return ((*calcEnergyFun)(x_ij)).cast<double>();
}

void PotentialOrder2Wrapper::calculateForce(model::Vec3 &force, const model::Vec3 &x_ij) const {
    pybind11::gil_scoped_acquire lock;
    force += ((*calcForceFun)(x_ij)).cast<readdy::model::Vec3>();
}

void
PotentialOrder2Wrapper::calculateForceAndEnergy(model::Vec3 &force, double &energy, const model::Vec3 &x_ij) const {
    pybind11::gil_scoped_acquire lock;
    energy += ((*calcEnergyFun)(x_ij)).cast<double>();
    force += ((*calcForceFun)(x_ij)).cast<readdy::model::Vec3>();
}

double PotentialOrder2Wrapper::getCutoffRadiusSquared() const {
    return getCutoffRadius() * getCutoffRadius();
}


}


}

