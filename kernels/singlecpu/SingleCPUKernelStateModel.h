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
#include <readdy/model/Vec3.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {

            struct ParticleData {
                std::shared_ptr<std::vector<boost::uuids::uuid>> ids;
                std::shared_ptr<std::vector<readdy::model::Vec3>> positions;
                std::unique_ptr<std::vector<unsigned int>> type;

                ParticleData();
                void addParticles(const std::vector<model::Particle> particles);
            };

            class SingleCPUKernelStateModel : public readdy::model::KernelStateModel {
            public:
                virtual const model::time_step_type getCurrentTimeStep() const override;

                virtual void updateModel(readdy::model::time_step_type t, bool forces, bool distances) override;

                virtual void addParticle(const model::Particle &p) override;
                virtual void addParticles(const std::vector<model::Particle> &p) override;

                virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;


                SingleCPUKernelStateModel();
                ~SingleCPUKernelStateModel();
                // move
                SingleCPUKernelStateModel(SingleCPUKernelStateModel &&rhs);
                SingleCPUKernelStateModel& operator=(SingleCPUKernelStateModel &&rhs);

                ParticleData* getParticleData() const;

            private:
                struct Impl;
                std::unique_ptr<Impl> pimpl;
            };

        }
    }
}

#endif //READDY2_MAIN_SINGLECPUKERNELSTATEMODEL_H
