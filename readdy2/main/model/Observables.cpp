/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#include <readdy/model/Observables.h>

void readdy::model::ParticlePositionObservable::evaluate(readdy::model::time_step_type t) {
    if(t == t_current) return;
    t_current = t;
    std::vector<Vec3> result = kernel->getKernelStateModel().getParticlePositions();
    *ParticlePositionObservable::result = result;
}

