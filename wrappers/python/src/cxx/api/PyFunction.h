/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
