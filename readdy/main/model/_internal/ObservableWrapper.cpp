#include <readdy/common/Types.h>
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


void readdy::model::ObservableWrapper::operator()(readdy::model::time_step_type t) {
    callback(t);
}

readdy::model::ObservableWrapper::ObservableWrapper(readdy::model::Kernel *const kernel,
                                                    const readdy::model::ObservableType &observable,
                                                    unsigned int stride)
        : ObservableBase(kernel, stride), observable(observable) {
}

void readdy::model::ObservableWrapper::evaluate() {
    observable(t_current);
}









