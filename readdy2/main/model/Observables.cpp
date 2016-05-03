/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#include <readdy/model/Observables.h>
#include <readdy/model/Kernel.h>

namespace readdy {
    namespace model {
        void ParticlePositionObservable::evaluate(time_step_type t) {
            if(t == t_current) return;
            t_current = t;
            std::vector<Vec3> result = kernel->getKernelStateModel().getParticlePositions();
            *ParticlePositionObservable::result = result;
        }
        
    }
}

