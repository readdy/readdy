/**
 * << detailed description >>
 *
 * @file Programs.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */
#include <boost/python.hpp>
#include <readdy/kernel/singlecpu/programs/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <PyConverters.h>

namespace bpy = boost::python;

using core_add_particle = readdy::model::programs::AddParticle;

using prog_factory_t = readdy::model::programs::ProgramFactory;

using scpu_add_particle = readdy::kernel::singlecpu::programs::SingleCPUAddParticleProgram;


void exportPrograms() {
    auto fun = &prog_factory_t::createProgram<scpu_add_particle>;

    bpy::class_<prog_factory_t>("ProgramFactory", bpy::no_init)
            .def("create_program_add_particles", readdy::py::adapt_unique(fun));

    bpy::class_<scpu_add_particle, std::unique_ptr<scpu_add_particle>>("AddParticle", bpy::init<readdy::kernel::singlecpu::SingleCPUKernel *>())
            .def("execute", &scpu_add_particle::execute);
}