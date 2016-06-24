/**
 * << detailed description >>
 *
 * @file PyFunction.h
 * @brief << brief description >>
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
        template<typename R, typename... Args>
        struct PyFunction {
            PyFunction(const boost::python::object& object) : py_obj(new boost::python::object(object), [](boost::python::object* o) { interpreter_lock lock; delete o; }){
                fun = [&] (Args&&... args) -> R {
                    interpreter_lock lock;
                    return (*py_obj)(std::forward<Args>(args)...);
                };
            }

            R operator()(Args&&... args) {
                return fun(std::forward<Args>(args)...);
            }

        protected:
            std::function<R(Args...)> fun;
            std::shared_ptr<boost::python::object> py_obj;
        };

        template<typename... Args>
        struct PyFunction<void, Args...> {
            PyFunction(const boost::python::object& object) : py_obj(new boost::python::object(object), [](boost::python::object* o) { interpreter_lock lock; delete o; }){
                fun = [&] (Args&&... args) {
                    interpreter_lock lock;
                    (*py_obj)(std::forward<Args>(args)...);
                };
            }

            void operator()(Args&&... args) {
                return fun(std::forward<Args>(args)...);
            }

        protected:
            std::function<void(Args...)> fun;
            std::shared_ptr<boost::python::object> py_obj;
        };
    }
}

#endif //READDY_MAIN_PYFUNCTION_H
