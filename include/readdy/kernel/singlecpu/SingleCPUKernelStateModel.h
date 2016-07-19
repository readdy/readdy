/**
 * << detailed description >>
 *
 * @file SingleCPUKernelStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_SINGLECPUKERNELSTATEMODEL_H
#define READDY_MAIN_SINGLECPUKERNELSTATEMODEL_H


#include <readdy/model/KernelStateModel.h>
#include <memory>
#include <readdy/model/Vec3.h>
#include <readdy/kernel/singlecpu/model/SingleCPUParticleData.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {

            class SingleCPUKernelStateModel : public readdy::model::KernelStateModel {
            public:


                virtual void updateNeighborList() override;
                virtual void calculateForces() override;

                virtual void addParticle(const readdy::model::Particle &p) override;
                virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;

                virtual void removeParticle(const readdy::model::Particle &p) override;
                virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;

                virtual double getEnergy() const override;

                void increaseEnergy(double increase);

                SingleCPUKernelStateModel(readdy::model::KernelContext const* context);
                ~SingleCPUKernelStateModel();
                // move
                SingleCPUKernelStateModel(SingleCPUKernelStateModel &&rhs);
                SingleCPUKernelStateModel& operator=(SingleCPUKernelStateModel &&rhs);

                readdy::kernel::singlecpu::model::SingleCPUParticleData* getParticleData() const;

                void setNeighborList(std::unique_ptr<model::SingleCPUNeighborList> ptr);
                const model::SingleCPUNeighborList * const getNeighborList() const;

                virtual const std::vector<readdy::model::Particle> getParticles() const override;

            private:
                struct Impl;
                std::unique_ptr<Impl> pimpl;
            };

        }
    }
}

#endif //READDY_MAIN_SINGLECPUKERNELSTATEMODEL_H
