/**
 * << detailed description >>
 *
 * @file CPUStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <readdy/model/KernelStateModel.h>

#ifndef READDY_MAIN_CPUSTATEMODEL_H
#define READDY_MAIN_CPUSTATEMODEL_H

namespace readdy {
    namespace kernel {
        namespace cpu {
            class CPUStateModel : public readdy::model::KernelStateModel {


            public:
                virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;
                virtual const std::vector<readdy::model::Particle> getParticles() const override;
                virtual void updateNeighborList() override;
                virtual void calculateForces() override;
                virtual void addParticle(const readdy::model::Particle &p) override;
                virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;
                virtual void removeParticle(const readdy::model::Particle &p) override;
                virtual double getEnergy() const override;
            };
        }
    }
}

#endif //READDY_MAIN_CPUSTATEMODEL_H
