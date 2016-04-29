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

            ParticlePositionObservable(Kernel *const kernel, unsigned int stride = 1) : Observable(kernel, stride) { }

            virtual ~ParticlePositionObservable() {
            }

            virtual void evaluate() override;
        };
    }
}

#endif //READDY2_MAIN_OBSERVABLES_H
