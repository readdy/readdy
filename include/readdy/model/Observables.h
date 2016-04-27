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

#include <readdy/model/Observable.h>
#include <readdy/model/Vec3.h>
#include <vector>

namespace readdy {
    namespace model {
        class ParticlePositionObservable : public Observable {
        public:
            ParticlePositionObservable(unsigned int stride, signal_t &signal) : Observable(stride, signal) {
            }

            virtual ~ParticlePositionObservable() {
            }

        protected:
            virtual void _evaluate(const std::shared_ptr<KernelContext> &context, const std::shared_ptr<KernelStateModel> &model) override;
        };
    }
}

#endif //READDY2_MAIN_OBSERVABLES_H
