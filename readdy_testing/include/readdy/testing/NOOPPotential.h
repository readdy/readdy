//
// Created by chris on 18.08.16.
//

#ifndef READDY_MAIN_NOOPPOTENTIAL_H
#define READDY_MAIN_NOOPPOTENTIAL_H

#include <readdy/model/potentials/PotentialOrder2.h>

namespace readdy {
    namespace testing {

        struct NOOPPotential : public readdy::model::potentials::PotentialOrder2 {
            NOOPPotential() : PotentialOrder2("no op") {}

            virtual double getCutoffRadius() const override {
                return 5;
            }

            virtual double calculateEnergy(const readdy::model::Vec3 &x_ij) const override { return 20; }

            virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) const override {}

            virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &x_ij) const override {}

            virtual double getMaximalForce(double kbt) const noexcept override { return 42; }

            virtual NOOPPotential *replicate() const override { return new NOOPPotential(*this); }
        };
    }
}

#endif //READDY_MAIN_NOOPPOTENTIAL_H
