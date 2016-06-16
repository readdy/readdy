/**
 * << detailed description >>
 *
 * @file interpreter_lock.h
 * @brief << brief description >>
 * @author clonker
 * @date 27.04.16
 */

#ifndef READDY_MAIN_INTERPRETER_LOCK_H
#define READDY_MAIN_INTERPRETER_LOCK_H

#include <Python.h>
#include <pystate.h>

namespace readdy {
    namespace py {
        class interpreter_lock {
        public:
            interpreter_lock();
            ~interpreter_lock();
        private:
            PyGILState_STATE gilState;
        };
    }
}


#endif //READDY_MAIN_INTERPRETER_LOCK_H
