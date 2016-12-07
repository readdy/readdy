/**
 * << detailed description >>
 *
 * @file TrajectoryObservable.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.08.16
 */

#ifndef READDY_MAIN_TRAJECTORYOBSERVABLE_H
#define READDY_MAIN_TRAJECTORYOBSERVABLE_H

#include <readdy/model/observables/Observables.h>

namespace readdy {
    namespace io {
        class TrajectoryWriter : public readdy::model::observables::Combiner<void, readdy::model::observables::ParticlePosition> {
            using kernel_t = readdy::model::Kernel;
            using ppObs_t = readdy::model::observables::ParticlePosition;
        public:
            TrajectoryWriter(const std::string& path, kernel_t *const kernel, unsigned int stride, ppObs_t* parent)
                    : readdy::model::observables::Combiner(kernel, stride, parent) {};

            virtual void evaluate() = 0;

        protected:

        };
    }
}
#endif //READDY_MAIN_TRAJECTORYOBSERVABLE_H
