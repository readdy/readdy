/**
 * << detailed description >>
 *
 * @file CPUDefaultReactionProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/UncontrolledApproximation.h>

using particle_t = readdy::model::Particle;

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                namespace reactions {
                    UncontrolledApproximation::UncontrolledApproximation(const CPUKernel *const kernel)
                            : kernel(kernel) {

                    }

                    void UncontrolledApproximation::execute() {
                    }
                }

            }
        }
    }
}