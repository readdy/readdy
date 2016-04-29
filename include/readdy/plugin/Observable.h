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

namespace readdy {
    namespace plugin {

        class Observable {
        public:

            Observable(unsigned int stride) : stride(stride) {
            };

            void setStride(const unsigned int stride) {
                Observable::stride = stride;
            }

            Observable() : Observable(1) {};

            virtual ~Observable() {
            };

            virtual void evaluate(const std::shared_ptr<readdy::model::KernelContext> &context, const std::shared_ptr<readdy::model::KernelStateModel> &model) = 0;

        protected:
            unsigned int stride;
        };
    }
}
#endif //READDY2_MAIN_OBSERVABLE_H
