/**
 * << detailed description >>
 *
 * @file Observables.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#ifndef READDY2_MAIN_OBSERVABLES_H
#define READDY2_MAIN_OBSERVABLES_H

#include <readdy/plugin/Observable.h>
#include <readdy/model/Vec3.h>
#include <vector>

namespace readdy {
    namespace plugin {
        class ParticlePositionObservable : public Observable {
        public:
            ParticlePositionObservable(unsigned int stride) : Observable(stride) {
            }

            virtual ~ParticlePositionObservable() {
            }

            virtual void evaluate(const std::shared_ptr<readdy::model::KernelContext> &context, const std::shared_ptr<readdy::model::KernelStateModel> &model) override;
        };
    }
}

#endif //READDY2_MAIN_OBSERVABLES_H
