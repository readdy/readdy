/**
 * SingleCPU declarations of potentials of second order (particle-particle interactions).
 *
 * @file PotentialsOrder2.h
 * @brief Declare potentials of second order.
 * @author clonker
 * @date 09.06.16
 */

#ifndef READDY_MAIN_POTENTIALSORDER2_H_H
#define READDY_MAIN_POTENTIALSORDER2_H_H

#include <readdy/model/potentials/PotentialsOrder2.h>
#include <readdy/model/Kernel.h>

namespace readdy {

namespace kernel {
namespace scpu {

namespace potentials {
class HarmonicRepulsion : public readdy::model::potentials::HarmonicRepulsion {
    using vec_t = readdy::model::Vec3;

public:
    HarmonicRepulsion(const readdy::model::Kernel *const kernel);

    virtual double calculateEnergy(const vec_t &x_ij) const override;

    virtual void calculateForce(vec_t &force, const vec_t &x_ij) const override;

    virtual void calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_ij) const override;

    virtual double getCutoffRadius() const override;

    virtual double getCutoffRadiusSquared() const override;

};

class WeakInteractionPiecewiseHarmonic : public readdy::model::potentials::WeakInteractionPiecewiseHarmonic {
    using vec_t = readdy::model::Vec3;
public:
    WeakInteractionPiecewiseHarmonic(const readdy::model::Kernel *const kernel);

    virtual double calculateEnergy(const readdy::model::Vec3 &x_ij) const override;

    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) const override;

    virtual void
    calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &x_ij) const override;

    virtual double getCutoffRadius() const override;

    virtual double getCutoffRadiusSquared() const override;

};

}
}
}
}
#endif //READDY_MAIN_POTENTIALSORDER2_H_H
