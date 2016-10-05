/**
 * << detailed description >>
 *
 * @file py_observable.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.06.16
 */

#ifndef READDY_MAIN_PY_OBSERVABLE_H
#define READDY_MAIN_PY_OBSERVABLE_H

#include <pybind11/pybind11.h>

#include <readdy/model/Kernel.h>

namespace readdy {
namespace py {
class PyObservable : public readdy::model::ObservableBase {

public:
    PyObservable(readdy::model::Kernel *const kernel, unsigned int stride, const pybind11::object &observableFun);

    virtual void evaluate() override;

private:
    std::shared_ptr<pybind11::object> py_ptr;
};
}
}


#endif //READDY_MAIN_PY_OBSERVABLE_H
