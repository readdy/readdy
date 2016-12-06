/**
 * << detailed description >>
 *
 * @file P1Cube.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder1.h>

namespace readdy {

namespace kernel {
namespace scpu {
namespace potentials {
SCPUCubePotential::SCPUCubePotential(const readdy::model::Kernel *const kernel) : readdy::model::potentials::CubePotential(
        kernel) {}

double SCPUCubePotential::calculateEnergy(const readdy::model::Vec3 &position) const {

    auto r = particleRadius;
    if (!isConsiderParticleRadius()) r = 0;

    double energy = 0;

    for (auto i = 0; i < 3; ++i) {
        if (position[i] - r < min[i] || position[i] + r > max[i]) {
            if (position[i] - r < min[i]) {
                energy += 0.5 * forceConstant * (position[i] - r - min[i]) * (position[i] - r - min[i]);
            } else {
                energy += 0.5 * forceConstant * (position[i] + r - max[i]) * (position[i] + r - max[i]);
            }
        }
    }

    return energy;
}

void SCPUCubePotential::calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) const {

    auto r = particleRadius;
    if (!isConsiderParticleRadius()) r = 0;
    for (auto i = 0; i < 3; i++) {
        if (position[i] - r < min[i] || position[i] + r > max[i]) {
            if (position[i] - r < min[i]) {
                force[i] += -1 * forceConstant * (position[i] - r - min[i]);
            } else {
                force[i] += -1 * forceConstant * (position[i] + r - max[i]);
            }
        }
    }

}

void SCPUCubePotential::calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                            const readdy::model::Vec3 &position) const {

    energy += calculateEnergy(position);
    calculateForce(force, position);

}

SCPUSpherePotential::SCPUSpherePotential(const readdy::model::Kernel *const kernel) : readdy::model::potentials::SpherePotential(kernel) { }

double SCPUSpherePotential::calculateEnergy(const readdy::model::Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    double energy = 0;
    if (distanceFromSphere > 0) {
        energy = 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
    }
    return energy;
}

void SCPUSpherePotential::calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere > 0) {
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

void SCPUSpherePotential::calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere > 0) {
        energy += 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}


}
}
}

}
