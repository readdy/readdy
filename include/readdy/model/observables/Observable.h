/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


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

#include <readdy/common/make_unique.h>
#include <readdy/common/signals.h>
#include <readdy/common/logging.h>
#include <readdy/common/Utils.h>

namespace readdy {
namespace model {
class Kernel;

namespace observables {
using time_step_type = unsigned long;
using signal_type = readdy::signals::signal<void(time_step_type)>;
using observable_type = signal_type::slot_type;

template<typename T>
struct ObservableName {
};


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

    const observables::time_step_type &getCurrentTimeStep() const {
        return t_current;
    }

    virtual ~ObservableBase() = default;

    virtual void callback(observables::time_step_type t) {
        if (should_execute_callback(t)) {
            firstCall = false;
            t_current = t;
            evaluate();
        }
    };

    virtual bool should_execute_callback(observables::time_step_type t) const {
        return (t_current != t || firstCall) && (stride == 0 || t % stride == 0);
    }

    virtual void evaluate() = 0;

protected:
    unsigned int stride;
    readdy::model::Kernel *const kernel;
    observables::time_step_type t_current = 0;
    bool firstCall = true;
};

template<typename Result>
class Observable : public ObservableBase {
public:
    using callback_function = std::function<void(const Result&)>;
    typedef Result result_t;

    Observable(Kernel *const kernel, unsigned int stride) : ObservableBase(kernel, stride), result() {
    }

    const result_t &getResult() {
        return result;
    }

    void setCallback(const callback_function& callbackFun) {
        Observable::_callback_f = std::move(callbackFun);
    }

    virtual void callback(observables::time_step_type t) override {
        if (should_execute_callback(t)) {
            ObservableBase::callback(t);
            _callback_f(result);
        }
    }


protected:
    Result result;
    callback_function _callback_f = [](const Result) {};
};

template<typename Res_t, typename... ParentObs_t>
class Combiner : public Observable<Res_t> {
public:
    Combiner(Kernel *const kernel, unsigned int stride, ParentObs_t*... parents)
            : Observable<Res_t>(kernel, stride), parentObservables(std::forward<ParentObs_t*>(parents)...) {}

    virtual void callback(observables::time_step_type t) override {
        readdy::util::collections::for_each_in_tuple(parentObservables, CallbackFunctor(ObservableBase::t_current));
        ObservableBase::callback(t);
    }

protected:
    std::tuple<ParentObs_t*...> parentObservables;
private:
    struct CallbackFunctor {
        observables::time_step_type currentTimeStep;

        CallbackFunctor(observables::time_step_type currentTimeStep) : currentTimeStep(currentTimeStep) {}

        template<typename T>
        void operator()(T *const obs) {
            if(obs->getCurrentTimeStep() != currentTimeStep) obs->callback(currentTimeStep);
        }
    };
};


}
}
}
#endif //READDY_MAIN_OBSERVABLE_H
