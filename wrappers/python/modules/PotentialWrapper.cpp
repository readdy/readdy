#include <readdy/model/Vec3.h>

/**
 * << detailed description >>
 *
 * @file PotentialWrapper.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.06.16
 */

#include "PotentialWrapper.h"

namespace readdy {
    double py::PotentialOrder2Wrapper::calculateEnergy(const model::Vec3 &x_i, const model::Vec3 &x_j) {
        return calcEnergyFun(x_i, x_j);
    }

    void py::PotentialOrder2Wrapper::calculateForce(model::Vec3 &force, const model::Vec3 &x_i, const model::Vec3 &x_j) {
        force += calcForceFun(x_i, x_j);
    }

    void py::PotentialOrder2Wrapper::calculateForceAndEnergy(model::Vec3 &force, double &energy, const model::Vec3 &x_i, const model::Vec3 &x_j) {
        energy += calcEnergyFun(x_i, x_j);
        force += calcForceFun(x_i, x_j);
    }


}

