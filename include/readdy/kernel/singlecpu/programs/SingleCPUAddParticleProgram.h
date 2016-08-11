/**
 * << detailed description >>
 *
 * @file SingleCPUAddParticleProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_SINGLECPUADDPARTICLEPROGRAM_H
#define READDY_MAIN_SINGLECPUADDPARTICLEPROGRAM_H

#include <readdy/model/programs/Programs.h>
#include <readdy/model/Particle.h>
#include <vector>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                class SingleCPUAddParticleProgram : public readdy::model::programs::AddParticle {
                public:
                    SingleCPUAddParticleProgram(SingleCPUKernel* kernel);

                    virtual ~SingleCPUAddParticleProgram() override;

                    virtual void execute() override;

                    void addParticle(readdy::model::Particle particle) {
                        particles.emplace_back(std::move(particle));
                    }

                    void setParticles(std::vector<readdy::model::Particle> particles) {
                        this->particles = std::move(particles) ;
                    }

                private:
                    std::vector<readdy::model::Particle> particles;
                    SingleCPUKernel *kernel;
                };
            }
        }
    }
}


#endif //READDY_MAIN_SINGLECPUADDPARTICLEPROGRAM_H
