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
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace readdy {
namespace rpy {
class PyProgram : public readdy::model::programs::Program {
    using super = readdy::model::programs::Program;
public:
    using super::Program;

    virtual void execute() override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE(void, readdy::model::programs::Program, execute,)
    }
};
}
}

#endif //READDY_MAIN_PROGRAMWRAP_H
