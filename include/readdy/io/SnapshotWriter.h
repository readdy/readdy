/**
 * << detailed description >>
 *
 * @file SnapshotObservable.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.08.16
 */

#ifndef READDY_MAIN_SNAPSHOTOBSERVABLE_H
#define READDY_MAIN_SNAPSHOTOBSERVABLE_H

#include <readdy/model/observables/Observable.h>

namespace readdy {
    namespace io {
        class SnapshotWriter : public readdy::model::observables::Observable<void> {
        public:
            SnapshotWriter(std::string fName,
                               model::Kernel *const kernel, unsigned int stride) : Observable(kernel, stride),
                                                                                   fileName(std::move(fName)) {
                // todo write config (ie kernel context)
            }
            virtual void evaluate() override {
                // todo write positions (virtual)
                // todo write forces (virtual)
            }
            
        protected:
            std::string fileName;
        };
    }
}

#endif //READDY_MAIN_SNAPSHOTOBSERVABLE_H
