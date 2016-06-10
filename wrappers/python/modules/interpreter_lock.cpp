/**
 * << detailed description >>
 *
 * @file interpreter_lock.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.04.16
 */

#include "interpreter_lock.h"

readdy::py::interpreter_lock::interpreter_lock() {
    gilState = PyGILState_Ensure();
}

readdy::py::interpreter_lock::~interpreter_lock() {
    PyGILState_Release(gilState);
}