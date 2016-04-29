/**
 * << detailed description >>
 *
 * @file Observable.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.04.16
 */

#ifndef READDY2_MAIN_OBSERVABLE_H
#define READDY2_MAIN_OBSERVABLE_H

#include <boost/signals2/signal.hpp>
#include <readdy/model/KernelContext.h>
#include <readdy/model/KernelStateModel.h>
#include "Kernel.h"

namespace readdy {
    namespace plugin {

        class Observable {
        public:

            Observable(Kernel *const kernel, unsigned int stride = 1) : stride(stride), kernel(kernel) {
            };

            void setStride(const unsigned int stride) {
                Observable::stride = stride;
            }

            virtual ~Observable() {
            };

            virtual void evaluate() = 0;

        protected:
            unsigned int stride;
            Kernel *const kernel;
        };
    }
}
#endif //READDY2_MAIN_OBSERVABLE_H
