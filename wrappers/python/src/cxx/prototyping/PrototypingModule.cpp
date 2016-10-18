/**
 * << detailed description >>
 *
 * @file PrototypingModule.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include "KernelWrap.h"
#include "../api/PyFunction.h"

namespace bpy = pybind11;

using rvp = bpy::return_value_policy;

void exportPrograms(bpy::module &);

void exportModelClasses(bpy::module &);

void exportPotentials(bpy::module &);

namespace scpu = readdy::kernel::singlecpu;

using rdy_scpu_model_t = scpu::SingleCPUKernelStateModel;
using scpu_kernel_t = scpu::SingleCPUKernel;

using scpu_kernel_wrap_t = scpu_kernel_t; // todo: do i need readdy::py::SingleCPUKernelWrap here?

using core_kernel_t = readdy::model::Kernel;
using core_kernel_wrap_t = readdy::py::KernelWrap;
using core_program_factory = readdy::model::programs::ProgramFactory;
using core_program_t = readdy::model::programs::Program;


// module
PYBIND11_PLUGIN (prototyping) {

    if(!readdy::log::console()) {
        spdlog::set_sync_mode();
        auto console = spdlog::stdout_color_mt("console");
        console->set_level(spdlog::level::debug);
        console->set_pattern("[          ] [%Y-%m-%d %H:%M:%S] [%t] [%l] %v");
    }

    bpy::module proto("prototyping", "ReaDDy prototyping python module");

    exportPrograms(proto);
    exportModelClasses(proto);
    exportPotentials(proto);

    bpy::class_<scpu_kernel_wrap_t>(proto, "SingleCPUKernel")
            .def(bpy::init<>())
            .def("get_kernel_state_model", &scpu_kernel_wrap_t::getKernelStateModel, rvp::reference_internal)
            .def("get_kernel_context", &scpu_kernel_wrap_t::getKernelContext, rvp::reference_internal)
            .def("get_available_potentials", &scpu_kernel_wrap_t::getAvailablePotentials)
            .def("get_potential_factory", &scpu_kernel_wrap_t::getPotentialFactory, rvp::reference_internal)
            .def("get_reaction_factory", &scpu_kernel_wrap_t::getReactionFactory, rvp::reference_internal)
            .def("get_observable_factory", &scpu_kernel_wrap_t::getObservableFactory, rvp::reference_internal)
            .def("get_program_factory", &scpu_kernel_wrap_t::getProgramFactory, rvp::reference_internal);

    return proto.ptr();
}