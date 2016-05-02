/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#include <readdy/model/Observables.h>

void readdy::model::ParticlePositionObservable::evaluate() {
    kernel->getKernelStateModel();
}

