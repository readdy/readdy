/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#include <readdy/model/Observables.h>

void readdy::model::ParticlePositionObservable::_evaluate(const std::shared_ptr<KernelContext> &context, const std::shared_ptr<KernelStateModel> &model) {
    model->getParticlePositions();
    // TODO
}

