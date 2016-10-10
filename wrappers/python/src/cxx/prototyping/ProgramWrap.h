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

#include <readdy/model/programs/Program.h>

namespace readdy {
namespace py {
class PyProgram : public readdy::model::programs::Program {
    using super = readdy::model::programs::Program;
public:
    using super::Program;

    virtual void execute() override {
        PYBIND11_OVERLOAD_PURE(void, readdy::model::programs::Program, execute,)
    }
};
}
}

#endif //READDY_MAIN_PROGRAMWRAP_H
