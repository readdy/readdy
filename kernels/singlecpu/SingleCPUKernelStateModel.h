/**
 * << detailed description >>
 *
 * @file SingleCPUKernelStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY2_MAIN_SINGLECPUKERNELSTATEMODEL_H
#define READDY2_MAIN_SINGLECPUKERNELSTATEMODEL_H


#include <readdy/model/KernelStateModel.h>
#include <memory>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUKernelStateModel : public readdy::model::KernelStateModel {
            public:
                virtual void updateModel(bool forces, bool distances) override;

                virtual void addParticle(const model::Particle &p) override;
                virtual void addParticles(const std::vector<model::Particle> &p) override;

                virtual std::vector<std::array<double, 3>> getParticlePositions() override;


                std::shared_ptr<std::vector<readdy::model::Particle>> getParticles() const;

                SingleCPUKernelStateModel();
                ~SingleCPUKernelStateModel();
                // move
                SingleCPUKernelStateModel(SingleCPUKernelStateModel &&rhs);
                SingleCPUKernelStateModel& operator=(SingleCPUKernelStateModel &&rhs);

            private:
                struct Impl;
                std::unique_ptr<Impl> pimpl;
            };

        }
    }
}

#endif //READDY2_MAIN_SINGLECPUKERNELSTATEMODEL_H
