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

#include <boost/function.hpp>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <boost/python/object.hpp>

namespace readdy {
    namespace py {
        class PotentialOrder2Wrapper : public readdy::model::potentials::PotentialOrder2 {

        public:
            PotentialOrder2Wrapper(const std::string &name, boost::function<double(model::Vec3, model::Vec3)> calculateEnergy, boost::function<model::Vec3(model::Vec3, model::Vec3)> calculateForce)
                    : PotentialOrder2(name), calcEnergyFun(calculateEnergy), calcForceFun(calculateForce) { }

            virtual double calculateEnergy(const model::Vec3 &x_i, const model::Vec3 &x_j) override;
            virtual void calculateForce(model::Vec3 &force, const model::Vec3 &x_i, const model::Vec3 &x_j) override;
            virtual void calculateForceAndEnergy(model::Vec3 &force, double &energy, const model::Vec3 &x_i, const model::Vec3 &x_j) override;

        protected:
            boost::function<double(model::Vec3, model::Vec3)> calcEnergyFun;
            boost::function<model::Vec3(model::Vec3, model::Vec3)> calcForceFun;
        };
    }
}

#endif //READDY_MAIN_POTENTIALWRAPPER_H
