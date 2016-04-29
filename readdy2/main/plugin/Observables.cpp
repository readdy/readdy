/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#include <readdy/plugin/Observables.h>

void readdy::plugin::ParticlePositionObservable::evaluate() {
    kernel->getKernelStateModel();
}

