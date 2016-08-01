/**
 * << detailed description >>
 *
 * @file PotentialWrapper.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.06.16
 */

#include "PyPotential.h"
#include "interpreter_lock.h"
#include <readdy/model/Vec3.h>
#include <iostream>
#include <boost/python/extract.hpp>

namespace readdy {

    namespace py {
        PotentialOrder2Wrapper::PotentialOrder2Wrapper(const std::string &name, boost::python::object o1, boost::python::object o2)
                : PotentialOrder2(name),
                  calcEnergyFun(new boost::python::object(o1), [](boost::python::object *o) {
                      interpreter_lock lock;
                      delete o;
                  }),
                  calcForceFun(new boost::python::object(o2), [](boost::python::object *o) {
                      interpreter_lock lock;
                      delete o;
                  })
                  { };

        double PotentialOrder2Wrapper::calculateEnergy(const model::Vec3 &x_ij) const {
            interpreter_lock lock;
            return boost::python::extract<double>((*calcEnergyFun)(x_ij));
        }

        void PotentialOrder2Wrapper::calculateForce(model::Vec3 &force, const model::Vec3 &x_ij) const {
            interpreter_lock lock;
            force += boost::python::extract<readdy::model::Vec3>((*calcForceFun)(x_ij));
        }

        void PotentialOrder2Wrapper::calculateForceAndEnergy(model::Vec3 &force, double &energy, const model::Vec3 &x_ij) const {
            interpreter_lock lock;
            energy += boost::python::extract<double>((*calcEnergyFun)(x_ij));
            force += boost::python::extract<readdy::model::Vec3>((*calcForceFun)(x_ij));
        }

        PotentialOrder2Wrapper *PotentialOrder2Wrapper::replicate() const {
            return new PotentialOrder2Wrapper(*this);
        }


    }


}

