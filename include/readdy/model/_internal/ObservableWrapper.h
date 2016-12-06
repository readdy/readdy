/**
 * This header file contains the definition of ObservableWrapper. It wraps a function object of the form
 * std::function<void(readdy::model::time_step_type)> together with a pointer to the kernel and a stride attribute
 * into an observable.
 *
 * @file ObservableWrapper.h
 * @brief This header file contains the definition of ObservableWrapper, which wraps a std::function to an observable.
 * @author clonker
 * @date 10.05.16
 */

#ifndef READDY_MAIN_OBSERVABLEWRAPPER_H
#define READDY_MAIN_OBSERVABLEWRAPPER_H

#include <readdy/model/observables/Observable.h>

namespace readdy {
namespace model {
namespace observables {

class ObservableWrapper : public ObservableBase {
public:
    ObservableWrapper(Kernel *const kernel,
                      const observables::observable_type &observable, unsigned int stride = 1);

    void operator()(observables::time_step_type t);

    virtual void evaluate() override;

protected:
    const observables::observable_type observable;
};

}
}
}
#endif //READDY_MAIN_OBSERVABLEWRAPPER_H
