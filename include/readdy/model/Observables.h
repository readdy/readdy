/**
 * << detailed description >>
 *
 * @file Observables.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#ifndef READDY_MAIN_OBSERVABLES_H
#define READDY_MAIN_OBSERVABLES_H

#include <readdy/model/Observable.h>
#include <readdy/model/Vec3.h>
#include <vector>
#include <readdy/common/Types.h>

namespace readdy {
    namespace model {


        class ParticlePositionObservable : public Observable<std::vector<Vec3>> {
        public:


            ParticlePositionObservable(Kernel *const kernel, unsigned int stride = 1) : Observable(kernel, stride) { }

            virtual ~ParticlePositionObservable() {
            }

            virtual void evaluate() override;
        };
        template<> struct ObservableName<ParticlePositionObservable> { static const std::string value; };

        class TestCombinerObservable : public CombinerObservable<std::vector<double>, ParticlePositionObservable, ParticlePositionObservable> {
        public:

            TestCombinerObservable(Kernel *const kernel, ParticlePositionObservable * obs1, ParticlePositionObservable * obs2, unsigned int stride)
                    : CombinerObservable(kernel, obs1, obs2, stride) {
            }

            virtual void evaluate() override;


        };
    }
}

#endif //READDY_MAIN_OBSERVABLES_H
