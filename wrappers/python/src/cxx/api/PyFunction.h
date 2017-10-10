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

namespace readdy {
namespace rpy {

template<typename Signature>
struct PyFunction;

template<typename R, typename... Args>
struct PyFunction<R(Args...)> {
    PyFunction(pybind11::object object) : py_obj(new pybind11::object(object), [](pybind11::object *o) {
        pybind11::gil_scoped_acquire lock;
        delete o;
    }) {
    }

    PyFunction(const PyFunction&) = default;
    PyFunction& operator=(const PyFunction&) = default;

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
    explicit PyFunction(pybind11::object object) : py_obj(new pybind11::object(object), [](pybind11::object *o) {
        pybind11::gil_scoped_acquire lock;
        delete o;
    }) {
    }

    PyFunction(const PyFunction&) = default;
    PyFunction& operator=(const PyFunction&) = default;

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
