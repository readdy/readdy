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

namespace readdy {
    namespace kernel {
        namespace singlecpu {

            class SingleCPUKernelStateModel : public readdy::model::KernelStateModel {
            public:
                virtual const readdy::model::time_step_type getCurrentTimeStep() const override;

                virtual void updateModel(readdy::model::time_step_type t, bool forces, bool distances) override;

                virtual void addParticle(const readdy::model::Particle &p) override;
                virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;

                virtual void removeParticle(const readdy::model::Particle &p) override;
                virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;


                SingleCPUKernelStateModel();
                ~SingleCPUKernelStateModel();
                // move
                SingleCPUKernelStateModel(SingleCPUKernelStateModel &&rhs);
                SingleCPUKernelStateModel& operator=(SingleCPUKernelStateModel &&rhs);

                readdy::kernel::singlecpu::model::SingleCPUParticleData* getParticleData() const;

            private:
                struct Impl;
                std::unique_ptr<Impl> pimpl;
            };

        }
    }
}

#endif //READDY_MAIN_SINGLECPUKERNELSTATEMODEL_H
