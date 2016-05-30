/**
 * << detailed description >>
 *
 * @file gil_lock.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.04.16
 */

#include "gil_lock.h"

readdy::py::gil_lock() {
    gilState = PyGILState_Ensure();
}

readdy::py::~gil_lock() {
    PyGILState_Release(gilState);
}