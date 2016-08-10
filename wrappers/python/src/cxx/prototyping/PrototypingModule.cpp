/**
 * << detailed description >>
 *
 * @file PrototypingModule.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <boost/python.hpp>
// singlecpu includes
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include "KernelWrap.h"
#include "../PyFunction.h"


#if PY_MAJOR_VERSION >= 3
int
#else

void
#endif
init_numpy() {
    if (PyArray_API == NULL) {
        import_array();
    }
}

void exportPrograms();
void exportModelClasses();

namespace bpy = boost::python;
namespace scpu = readdy::kernel::singlecpu;

using _rdy_scpu_model_t = scpu::SingleCPUKernelStateModel;
using scpu_kernel_t = scpu::SingleCPUKernel;

using scpu_kernel_wrap_t = scpu_kernel_t; // todo: do i need readdy::py::SingleCPUKernelWrap here?

using core_kernel_t = readdy::model::Kernel;
using core_kernel_wrap_t = readdy::py::KernelWrap;
using core_program_factory = readdy::model::programs::ProgramFactory;
using core_program_t = readdy::model::programs::Program;


// module
BOOST_PYTHON_MODULE (prototyping) {
    init_numpy();
    PyEval_InitThreads();

    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    bpy::docstring_options doc_options;
    doc_options.enable_all();

    exportPrograms();
    exportModelClasses();

    bpy::class_<scpu_kernel_wrap_t, boost::noncopyable>("SingleCPUKernel")
            .def("get_kernel_state_model", &scpu_kernel_wrap_t::getKernelStateModel, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("get_kernel_context", &scpu_kernel_wrap_t::getKernelContext, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("get_available_potentials", &scpu_kernel_wrap_t::getAvailablePotentials)
            .def("get_potential_factory", &scpu_kernel_wrap_t::getPotentialFactory, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("get_reaction_factory", &scpu_kernel_wrap_t::getReactionFactory, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("get_observable_factory", &scpu_kernel_wrap_t::getObservableFactory, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("get_program_factory", &scpu_kernel_wrap_t::getProgramFactory, bpy::return_value_policy<bpy::reference_existing_object>());
}