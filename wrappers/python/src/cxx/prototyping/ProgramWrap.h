/**
 * << detailed description >>
 *
 * @file ProgramWrap.h
 * @brief << brief description >>
 * @author clonker
 * @date 08.08.16
 */

#ifndef READDY_MAIN_PROGRAMWRAP_H
#define READDY_MAIN_PROGRAMWRAP_H

#include <boost/python.hpp>
#include <readdy/model/programs/Program.h>

namespace bpy = boost::python;

namespace readdy {
    namespace py {
        struct Program : public readdy::model::programs::Program, bpy::wrapper<readdy::model::programs::Program> {
            virtual void execute() override {
                this->get_override("execute");
            }
        };
    }
}

#endif //READDY_MAIN_PROGRAMWRAP_H
