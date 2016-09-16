/**
 * Second order potential which is only useful for testing. The behavior is defined by the parameters cutoff, force and energy
 * which the constructor accepts. This approach is similar to mocking, but actually mocking a PotentialOrder2 has too many side effects.
 *
 * @file NOOPPotential.h
 * @brief Second-order potential that has no real effect. Just returns values defined during construction.
 * @author chrisfroe
 * @date 19.08.16
 */

#ifndef READDY_MAIN_NOOPPOTENTIAL_H
#define READDY_MAIN_NOOPPOTENTIAL_H

#include <readdy/model/potentials/PotentialOrder2.h>

namespace readdy {
namespace testing {
struct NOOPPotentialOrder2 : public readdy::model::potentials::PotentialOrder2 {
    NOOPPotentialOrder2(double cutoff = 0, double force = 0, double energy = 0) : PotentialOrder2("no op"),
                                                                                  cutoff(cutoff), force(force),
                                                                                  energy(energy) {}

    virtual double getCutoffRadius() const override {
        return cutoff;
    }

    virtual double calculateEnergy(const readdy::model::Vec3 &x_ij) const override {
        return energy;
    }

    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) const override {
        force[0] = NOOPPotentialOrder2::force;
        force[1] = NOOPPotentialOrder2::force;
        force[2] = NOOPPotentialOrder2::force;
    }

    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                         const readdy::model::Vec3 &x_ij) const override {
        energy = NOOPPotentialOrder2::calculateEnergy(x_ij);
        NOOPPotentialOrder2::calculateForce(force, x_ij);
    }

    virtual double getMaximalForce(double kbt) const noexcept override { return force; }

    virtual double getCutoffRadiusSquared() const override {
        return getCutoffRadius() * getCutoffRadius();
    }

    virtual NOOPPotentialOrder2 *replicate() const override { return new NOOPPotentialOrder2(*this); }

    double cutoff, force, energy;
};
}
}

#endif //READDY_MAIN_NOOPPOTENTIAL_H
