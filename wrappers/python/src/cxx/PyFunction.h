/**
 * PyFunction provides a way to insert python function pointer into C++ code. It is
 * templated with respect to the result of the callback and the arguments.
 *
 * @file PyFunction.h
 * @brief PyFunction is a wrapper for a Python callable object that can be used as std::function object.
 * @author clonker
 * @date 24.06.16
 */

#ifndef READDY_MAIN_PYFUNCTION_H
#define READDY_MAIN_PYFUNCTION_H

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <utility>

namespace readdy {
namespace py {

template<typename Signature>
struct PyFunction;

template<typename R, typename... Args>
struct PyFunction<R(Args...)> {
    PyFunction(pybind11::object object) : py_obj(new pybind11::object(object), [](pybind11::object *o) {
        pybind11::gil_scoped_acquire lock;
        delete o;
    }) {
    }

    R operator()(Args &&... args) {
        pybind11::gil_scoped_acquire lock;
        pybind11::object res = (*py_obj)(std::forward<Args>(args)...);
        return res.cast<R>();
    }

protected:
    std::shared_ptr<pybind11::object> py_obj;
};

template<typename... Args>
struct PyFunction<void(Args...)> {
    PyFunction(pybind11::object object) : py_obj(new pybind11::object(object), [](pybind11::object *o) {
        pybind11::gil_scoped_acquire lock;
        delete o;
    }) {
    }

    void operator()(Args &&... args) {
        pybind11::gil_scoped_acquire lock;
        (*py_obj)(std::forward<Args>(args)...);
    }

protected:
    std::shared_ptr<pybind11::object> py_obj;
};

}
}

#endif //READDY_MAIN_PYFUNCTION_H
