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

#include <pybind11/pybind11.h>

#include <functional>
#include <readdy/model/potentials/PotentialOrder2.h>

namespace readdy {
namespace py {
class PotentialOrder2Wrapper : public readdy::model::potentials::PotentialOrder2 {

public:
    PotentialOrder2Wrapper(const std::string &name, pybind11::object o1, pybind11::object o2);

    virtual double calculateEnergy(const model::Vec3 &x_ij) const override;

    virtual void calculateForce(model::Vec3 &force, const model::Vec3 &x_ij) const override;

    virtual void calculateForceAndEnergy(model::Vec3 &force, double &energy, const model::Vec3 &x_ij) const override;

    virtual double getCutoffRadius() const override {
        // todo!
        return 50;
    }

    virtual double getMaximalForce(double kbt) const noexcept override {
        // todo!
        return 0;
    }

    virtual double getCutoffRadiusSquared() const override;


protected:
    std::shared_ptr<pybind11::object> calcEnergyFun;
    std::shared_ptr<pybind11::object> calcForceFun;
};
}
}

#endif //READDY_MAIN_POTENTIALWRAPPER_H
