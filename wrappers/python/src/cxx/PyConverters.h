/**
 * << detailed description >>
 *
 * @file PyConverters.h
 * @brief << brief description >>
 * @author clonker
 * @date 28.06.16
 */

#ifndef READDY_MAIN_PYCONVERTERS_H
#define READDY_MAIN_PYCONVERTERS_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python.hpp>
#include <iostream>
#include <vector>
#include <list>

namespace readdy {
namespace py {

template<typename T, typename C, typename... Args, typename... Args2>
boost::python::object adapt_function(const std::function<T(Args2...)> &(C::*fn)(Args...) const) {
    return boost::python::make_function(
            [fn](C &self, Args... args) {
                return boost::python::make_function(
                        (self.*fn)(args...),
                        boost::python::return_value_policy<boost::python::reference_existing_object>(),
                        boost::mpl::vector<T, Args2...>()
                );
            }, boost::python::default_call_policies(), boost::mpl::vector<boost::python::object, C &>()
    );
};

template<typename T, typename C, typename... Args, typename... Args2>
boost::python::object adapt_function(const std::function<T(Args2...)> (*fn)(C, Args...)) {
    return boost::python::make_function(
            [fn](C &self, Args... args) {
                return boost::python::make_function(
                        (*fn)(self, args...),
                        boost::python::default_call_policies(),
                        boost::mpl::vector<T, Args2...>()
                );
            }, boost::python::default_call_policies(), boost::mpl::vector<boost::python::object, C &>()
    );
};

template<typename T, std::size_t N>
boost::python::list toList(std::array<T, N> arr) {
    boost::python::list list;
    for (T t : arr) {
        list.append<T>(t);
    }
    return list;
};

template<typename T>
std::list<T> sequence_to_list(const boost::python::object &o) {
    boost::python::stl_input_iterator<T> begin(o), end;
    return std::list<T>(begin, end);
}

template<typename T>
std::vector<T> sequence_to_vector(const boost::python::object &o) {
    boost::python::stl_input_iterator<T> begin(o), end;
    return std::vector<T>(begin, end);
}

/**
 * Takes a unique pointer and transfers its ownership to python
 * @param fn the function
 * @return a python object holding the unique pointer's contents
 */
template<typename T, typename C, typename ...Args>
boost::python::object adapt_unique(std::unique_ptr<T> (C::*fn)(Args...)) {
    return boost::python::make_function(
            [fn](C &self, Args... args) { return (self.*fn)(args...).release(); },
            boost::python::return_value_policy<boost::python::manage_new_object>(),
            boost::mpl::vector<T *, C &, Args...>()
    );
}

template<typename T, typename C, typename ...Args>
boost::python::object adapt_unique(std::unique_ptr<T> (C::*fn)(Args...) const) {
    return boost::python::make_function(
            [fn](C &self, Args... args) { return (self.*fn)(args...).release(); },
            boost::python::return_value_policy<boost::python::manage_new_object>(),
            boost::mpl::vector<T *, C &, Args...>()
    );
}

// Converts a std::pair instance to a Python tuple.
template<typename T1, typename T2>
struct std_pair_to_tuple {
    static PyObject *convert(std::pair<T1, T2> const &p) {
        return boost::python::incref(
                boost::python::make_tuple(p.first, p.second).ptr());
    }

    static PyTypeObject const *get_pytype() { return &PyTuple_Type; }
};

// Helper for convenience.
template<typename T1, typename T2>
struct std_pair_to_python_converter {
    std_pair_to_python_converter() {
        boost::python::to_python_converter<
                std::pair<T1, T2>,
                std_pair_to_tuple<T1, T2>,
                true //std_pair_to_tuple has get_pytype
        >();
    }
};

template<typename T>
struct std_vector_to_ndarray {
    static PyObject *convert(const std::vector<T> &vec) {
        npy_intp size = vec.size();
        double *data = size ? const_cast<double *>(&vec[0]) : static_cast<double *>(NULL);
        PyObject *obj = PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, data);
        boost::python::handle<> handle(obj);
        boost::python::numeric::array arr(handle);
        return boost::python::incref(arr.copy().ptr());
    }

    static PyTypeObject const *get_pytype() { return &PyArray_Type; }
};

template<typename T>
struct std_vector_to_python_converter {
    std_vector_to_python_converter() {
        boost::python::to_python_converter<
                std::vector<T>,
                std_vector_to_ndarray<T>,
                true //std_pair_to_tuple has get_pytype
        >();
    }
};
}
}


#endif //READDY_MAIN_PYCONVERTERS_H
