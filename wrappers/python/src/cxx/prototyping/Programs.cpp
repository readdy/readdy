/**
 * << detailed description >>
 *
 * @file Programs.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */
#include <boost/python.hpp>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>

namespace bpy = boost::python;

using add_particle = readdy::kernel::singlecpu::programs::SingleCPUAddParticleProgram;

void exportPrograms() {
    bpy::class_<add_particle>("AddParticle", bpy::init<readdy::kernel::singlecpu::SingleCPUKernel*>())
            .def("execute", &add_particle::execute)
            .def("set_particles", &add_particle::setParticles);
}