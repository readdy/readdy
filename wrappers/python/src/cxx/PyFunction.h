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

#include <utility>
#include <boost/python/object.hpp>
#include "interpreter_lock.h"

namespace readdy {
    namespace py {

        template<typename Signature>
        struct PyFunction;

        template<typename R, typename... Args>
        struct PyFunction<R(Args...)> {
            PyFunction(boost::python::object object) : py_obj(new boost::python::object(object), [](boost::python::object* o) { interpreter_lock lock; delete o; }){
            }

            R operator()(Args&&... args) {
                interpreter_lock lock;
                return boost::python::extract<R>((*py_obj)(std::forward<Args>(args)...));
            }

        protected:
            std::shared_ptr<boost::python::object> py_obj;
        };

        template<typename... Args>
        struct PyFunction<void(Args...)> {
            PyFunction(boost::python::object object) : py_obj(new boost::python::object(object), [](boost::python::object* o) { interpreter_lock lock; delete o; }){
            }

            void operator()(Args&&... args) {
                interpreter_lock lock;
                (*py_obj)(std::forward<Args>(args)...);
            }

        protected:
            std::shared_ptr<boost::python::object> py_obj;
        };

    }
}

#endif //READDY_MAIN_PYFUNCTION_H
