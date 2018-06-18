//
// Created by mho on 5/29/18.
//

#pragma once

#include <pybind11/numpy.h>

template<typename T>
using np_array = pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>;
