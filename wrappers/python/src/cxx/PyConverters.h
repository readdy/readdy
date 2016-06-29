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

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python.hpp>
#include <iostream>
#include <vector>

namespace readdy {
    namespace py {
        // Converts a std::pair instance to a Python tuple.
        template <typename T1, typename T2>
        struct std_pair_to_tuple
        {
            static PyObject* convert(std::pair<T1, T2> const& p)
            {
                return boost::python::incref(
                        boost::python::make_tuple(p.first, p.second).ptr());
            }
            static PyTypeObject const *get_pytype () {return &PyTuple_Type; }
        };

        // Helper for convenience.
        template <typename T1, typename T2>
        struct std_pair_to_python_converter
        {
            std_pair_to_python_converter()
            {
                boost::python::to_python_converter<
                        std::pair<T1, T2>,
                        std_pair_to_tuple<T1, T2>,
                        true //std_pair_to_tuple has get_pytype
                >();
            }
        };

        template <typename T>
        struct std_vector_to_ndarray {
            static PyObject* convert(const std::vector<T> &vec) {
                npy_intp size = vec.size();
                double * data = size ? const_cast<double *>(&vec[0]) : static_cast<double *>(NULL);
                PyObject* obj = PyArray_SimpleNewFromData( 1, &size, NPY_DOUBLE, data );
                boost::python::handle<> handle( obj );
                boost::python::numeric::array arr( handle );
                return boost::python::incref(arr.copy().ptr());
            }

            static PyTypeObject const *get_pytype () { return &PyArray_Type; }
        };

        template <typename T>
        struct std_vector_to_python_converter {
            std_vector_to_python_converter (){
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
