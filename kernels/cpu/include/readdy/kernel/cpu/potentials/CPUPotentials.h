/**
 * << detailed description >>
 *
 * @file CPUPotentials.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_MAIN_CPUPOTENTIALS_H
#define READDY_MAIN_CPUPOTENTIALS_H


#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/potentials/PotentialsOrder2.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace potentials {
                struct CPUCubePotential : public readdy::model::potentials::CubePotential {
                    CPUCubePotential(const readdy::model::Kernel *const kernel);
                    virtual double calculateEnergy(const readdy::model::Vec3 &position) override;
                    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) override;
                    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                                         const readdy::model::Vec3 &position) override;
                    virtual CPUCubePotential *replicate() const override;
                };

                struct CPUHarmonicRepulsion : public readdy::model::potentials::HarmonicRepulsion {
                    CPUHarmonicRepulsion(const readdy::model::Kernel *const kernel);
                    virtual double calculateEnergy(const readdy::model::Vec3 &x_ij) override;
                    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) override;
                    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                                         const readdy::model::Vec3 &x_ij) override;
                    virtual double getCutoffRadius() override;
                    virtual CPUHarmonicRepulsion *replicate() const override;
                };
            }
        }
    }
}
#endif //READDY_MAIN_CPUPOTENTIALS_H
