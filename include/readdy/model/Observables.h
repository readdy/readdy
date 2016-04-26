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

#include "Observable.h"
#include "Vec3.h"
#include <vector>

namespace readdy {
    namespace model {
        class ParticlePositionObservable : public Observable<std::vector<Vec3>()> {
        private:
            boost::signals2::connection connection;

            ParticlePositionObservable(int stride) : Observable(stride) {
                connection = signal.connect(std::bind(&readdy::model::ParticlePositionObservable::evaluate, this));
            }

            virtual ~ParticlePositionObservable() {
                connection.disconnect();
            }

        public:
            virtual std::vector<Vec3> evaluate() override;

        };
    }
}

#endif //READDY2_MAIN_OBSERVABLES_H
