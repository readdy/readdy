/**
 * << detailed description >>
 *
 * @file py_observable.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 14.06.16
 */

#include "PyObservable.h"
#include "interpreter_lock.h"

namespace readdy {
    namespace py {
        PyObservable::PyObservable(readdy::model::Kernel *const kernel, unsigned int stride, const boost::python::object &observableFun)
                : readdy::model::ObservableBase(kernel, stride),
                  py_ptr(new boost::python::object(observableFun),
                         [](boost::python::object *o) {
                             interpreter_lock lock;
                             delete o;
                         })
        {

        }

        void PyObservable::evaluate() {
            interpreter_lock lock;
            (*py_ptr)(*kernel);
        }
    }
}