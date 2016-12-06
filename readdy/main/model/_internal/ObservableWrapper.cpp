#include <readdy/model/_internal/ObservableWrapper.h>
#include <readdy/model/Kernel.h>

/**
 * << detailed description >>
 *
 * @file ObservableWrapper.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.05.16
 */


void readdy::model::observables::ObservableWrapper::operator()(observables::time_step_type t) {
    callback(t);
}

readdy::model::observables::ObservableWrapper::ObservableWrapper(readdy::model::Kernel *const kernel,
                                                                 const observables::observable_type &observable,
                                                                 unsigned int stride)
        : ObservableBase(kernel, stride), observable(observable) {
}

void readdy::model::observables::ObservableWrapper::evaluate() {
    observable(t_current);
}









