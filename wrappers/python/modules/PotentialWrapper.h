/**
 * << detailed description >>
 *
 * @file PotentialWrapper.h
 * @brief << brief description >>
 * @author clonker
 * @date 10.06.16
 */

#ifndef READDY_MAIN_POTENTIALWRAPPER_H
#define READDY_MAIN_POTENTIALWRAPPER_H

#include <functional>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <boost/python/object.hpp>

namespace readdy {
    namespace py {
        class PotentialOrder2Wrapper : public readdy::model::potentials::PotentialOrder2 {

        public:
            PotentialOrder2Wrapper(const std::string &name, boost::python::object o1, boost::python::object o2);

            virtual double calculateEnergy(const model::Vec3 &x_i, const model::Vec3 &x_j) override;
            virtual void calculateForce(model::Vec3 &force, const model::Vec3 &x_i, const model::Vec3 &x_j) override;
            virtual void calculateForceAndEnergy(model::Vec3 &force, double &energy, const model::Vec3 &x_i, const model::Vec3 &x_j) override;

        protected:
            std::shared_ptr<boost::python::object> calcEnergyFun;
            std::shared_ptr<boost::python::object> calcForceFun;
        };
    }
}

#endif //READDY_MAIN_POTENTIALWRAPPER_H
