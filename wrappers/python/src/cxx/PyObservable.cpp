/**
 * << detailed description >>
 *
 * @file py_observable.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 14.06.16
 */

#include "PyObservable.h"

namespace readdy {
namespace py {
PyObservable::PyObservable(readdy::model::Kernel *const kernel, unsigned int stride,
                           const pybind11::object &observableFun)
        : readdy::model::ObservableBase(kernel, stride),
          py_ptr(new pybind11::object(observableFun),
                 [](pybind11::object *o) {
                     pybind11::gil_scoped_acquire lock;
                     delete o;
                 }) {

}

void PyObservable::evaluate() {
    pybind11::gil_scoped_acquire lock;
    (*py_ptr)(*kernel);
}
}
}