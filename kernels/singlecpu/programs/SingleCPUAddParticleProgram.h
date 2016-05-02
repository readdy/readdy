/**
 * << detailed description >>
 *
 * @file SingleCPUAddParticleProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY2_MAIN_SINGLECPUADDPARTICLEPROGRAM_H
#define READDY2_MAIN_SINGLECPUADDPARTICLEPROGRAM_H

#include <readdy/plugin/Programs.h>
#include <readdy/model/Particle.h>
#include <vector>
#include "../SingleCPUKernel.h"

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                class SingleCPUAddParticleProgram : public readdy::plugin::AddParticleProgram {
                public:
                    SingleCPUAddParticleProgram(SingleCPUKernel* kernel);

                    virtual ~SingleCPUAddParticleProgram() override;

                    virtual void execute() override;

                    void addParticle(const readdy::model::Particle &particle) {
                        particles.emplace_back(particle);
                    }

                    void addParticle(readdy::model::Particle &&particle) {
                        particles.emplace_back(std::move(particle));
                    }

                    void setParticles(std::vector<readdy::model::Particle> &particles) {
                        this->particles = particles;
                    }

                private:
                    std::vector<readdy::model::Particle> particles;
                    SingleCPUKernel *kernel;
                };
            }
        }
    }
}


#endif //READDY2_MAIN_SINGLECPUADDPARTICLEPROGRAM_H