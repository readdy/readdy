/**
 * << detailed description >>
 *
 * @file CPUStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_MAIN_CPUSTATEMODEL_H
#define READDY_MAIN_CPUSTATEMODEL_H


#include <readdy/model/KernelStateModel.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/kernel/singlecpu/model/SingleCPUParticleData.h>
#include <readdy/kernel/cpu/model/NeighborList.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            class CPUStateModel : public readdy::model::KernelStateModel {

            public:
                CPUStateModel(readdy::model::KernelContext *const context);
                ~CPUStateModel();
                virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;
                virtual const std::vector<readdy::model::Particle> getParticles() const override;
                virtual void updateNeighborList() override;
                virtual void calculateForces() override;
                virtual void addParticle(const readdy::model::Particle &p) override;
                virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;
                virtual void removeParticle(const readdy::model::Particle &p) override;
                virtual double getEnergy() const override;

                readdy::kernel::singlecpu::model::SingleCPUParticleData *const getParticleData() const;

                const model::NeighborList*const getNeighborList() const;
            private:
                struct Impl;
                std::unique_ptr<Impl> pimpl;
            };
        }
    }
}

#endif //READDY_MAIN_CPUSTATEMODEL_H
