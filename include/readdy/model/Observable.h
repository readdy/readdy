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
#include <readdy/common/make_unique.h>
#include <readdy/common/Types.h>

#define DefineObservableName( observableName ) \
static std::string name()  \
{ \
     return #observableName; \
}

namespace readdy {
    namespace model {

        class Kernel;
        class Observable {
        public:

            Observable(readdy::model::Kernel *const kernel, unsigned int stride = 1) : stride(stride), kernel(kernel) {
            };

            void setStride(const unsigned int stride) {
                Observable::stride = stride;
            }

            unsigned int getStride() const {
                return stride;
            }

            virtual ~Observable() {
            };

            template<typename T>
            T getAs() {
                return dynamic_cast<T>(*this);
            }

            virtual void evaluate(readdy::model::time_step_type t) = 0;

        protected:
            unsigned int stride;
            readdy::model::Kernel *const kernel;
            readdy::model::time_step_type t_current;
        };

        template<typename Result>
        class ObservableWithResult : public Observable {
        public:
            ObservableWithResult(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride) {
                result = std::make_unique<Result>();
            }

            Result* getResult() {
                return result.get();
            }
        protected:
            std::unique_ptr<Result> result;
        };
    }
}
#endif //READDY2_MAIN_OBSERVABLE_H
