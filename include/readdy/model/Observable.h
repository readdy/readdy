/**
 * The struct ObservableName is specialized for every observable type, giving access to a static string that corresponds
 * to an observable's name.
 *
 * ObservableBase is the base class of all observables. They have a stride and a current time step. Also, they provide
 * an abstract evaluate method, in which the subclasses should execute their calculations.
 *
 * Observable derives from ObservableBase and is templateized to a certain result type.
 *
 * CombinerObservable takes two Observables and combines their results into a third result to avoid duplication of
 * work.
 *
 * @file Observable.h
 * @brief Header file containing the definitions for ObservableName, ObservableBase, Observable and CombinerObservable.
 * @author clonker
 * @date 22.04.16
 */

#ifndef READDY_MAIN_OBSERVABLE_H
#define READDY_MAIN_OBSERVABLE_H

#include <boost/signals2/signal.hpp>
#include <readdy/common/make_unique.h>
#include <readdy/common/Types.h>

namespace readdy {
    namespace model {
        template<typename T>
        struct ObservableName { };

        class Kernel;
        class ObservableBase {
        public:

            ObservableBase(readdy::model::Kernel *const kernel, unsigned int stride = 1) : stride(stride), kernel(kernel) {
            };

            void setStride(const unsigned int stride) {
                ObservableBase::stride = stride;
            }

            unsigned int getStride() const {
                return stride;
            }

            const readdy::model::time_step_type& getCurrentTimeStep() const {
                return t_current;
            }

            virtual ~ObservableBase();

            virtual void callback(readdy::model::time_step_type t) {
                if(t_current == t && !firstCall) return;
                firstCall = false;
                t_current = t;
                evaluate();
            };

            virtual void evaluate() = 0;

        protected:
            unsigned int stride;
            readdy::model::Kernel *const kernel;
            readdy::model::time_step_type t_current = 0;
            bool firstCall = true;
        };

        template<typename Result>
        class Observable : public ObservableBase {
        public:
            typedef Result result_t;
            Observable(Kernel *const kernel, unsigned int stride) : ObservableBase(kernel, stride) {
            }

            const result_t& getResult() {
                return result;
            }

            void setCallback(std::function<void(const result_t&)>&& callbackFun) {
                Observable::_callback_f = std::move(callbackFun);
            }

            void setCallback(const std::function<void(const result_t&)>& callbackFun) {
                Observable::_callback_f = callbackFun;
            }

            virtual void callback(readdy::model::time_step_type t) override {
                if(t_current == t && !firstCall) return;
                ObservableBase::callback(t);
                _callback_f(result);
            }


        protected:
            Result result;
            std::function<void(Result)> _callback_f = [](const Result){};
        };

        template<typename Res_t, typename Obs1_t, typename Obs2_t>
        class CombinerObservable : public Observable<Res_t> {
            static_assert(std::is_base_of<readdy::model::ObservableBase, Obs1_t>::value, "Type of Observable 1 was not a subtype of ObservableBase");
            static_assert(std::is_base_of<readdy::model::ObservableBase, Obs2_t>::value, "Type of Observable 2 was not a subtype of ObservableBase");
        public:
            typedef Obs1_t Observable1_type;
            typedef Obs2_t Observable2_type;
            CombinerObservable(Kernel *const kernel, Obs1_t * obs1, Obs2_t * obs2, unsigned int stride = 1)
                    : readdy::model::Observable<Res_t>::Observable(kernel, stride)
                    , obs1(obs1), obs2(obs2){ }

            virtual void callback(readdy::model::time_step_type t) override {
                if(obs1->getCurrentTimeStep() != ObservableBase::t_current) obs1->callback(ObservableBase::t_current);
                if(obs2->getCurrentTimeStep() != ObservableBase::t_current) obs2->callback(ObservableBase::t_current);
                ObservableBase::callback(t);
            }

        protected:
            Obs1_t * obs1;
            Obs2_t * obs2;
        };
    }
}
#endif //READDY_MAIN_OBSERVABLE_H
