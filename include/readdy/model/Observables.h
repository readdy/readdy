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
#include <readdy/common/Types.h>

namespace readdy {
    namespace model {

        template<typename Res_t, typename... Observables>
        class CombinerObservable : public Observable<Res_t> {
            CombinerObservable(Kernel *const kernel, unsigned int stride = 1) : readdy::model::Observable<Res_t>::Observable(kernel, stride) { }
        };

        class ParticlePositionObservable : public Observable<std::vector<Vec3>> {
        public:
            DefineObservableName(ParticlePositionObservable)

            ParticlePositionObservable(Kernel *const kernel, unsigned int stride = 1) : Observable(kernel, stride) { }

            virtual ~ParticlePositionObservable() {
            }

            virtual void evaluate(readdy::model::time_step_type t) override;
        };
    }
}

#endif //READDY2_MAIN_OBSERVABLES_H
