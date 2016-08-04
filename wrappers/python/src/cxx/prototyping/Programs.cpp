/**
 * << detailed description >>
 *
 * @file Programs.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */
#include <boost/python.hpp>
#include <PyConverters.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUUpdateNeighborList.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>

namespace bpy = boost::python;
namespace rp = readdy::py;

using prog_factory_t = readdy::model::programs::ProgramFactory;

using add_particle_t = readdy::kernel::singlecpu::programs::SingleCPUAddParticleProgram;
using euler_integrator_t = readdy::kernel::singlecpu::programs::SingleCPUEulerBDIntegrator;
using forces_t = readdy::kernel::singlecpu::programs::SingleCPUCalculateForces;
using neighbor_list_t = readdy::kernel::singlecpu::programs::SingleCPUUpdateNeighborList;

using reactions_uncontrolled_approx_t = readdy::kernel::singlecpu::programs::reactions::UncontrolledApproximation;

void exportPrograms() {
    auto f_add_particle = &prog_factory_t::createProgram<add_particle_t>;
    auto f_euler_integrator = &prog_factory_t::createProgram<euler_integrator_t>;
    auto f_forces = &prog_factory_t::createProgram<forces_t>;
    auto f_neighbor_list = &prog_factory_t::createProgram<neighbor_list_t>;
    auto f_reactions_uncontrolled_approximation = &prog_factory_t::createProgram<reactions_uncontrolled_approx_t>;

    using scpu_kernel_t = readdy::kernel::singlecpu::SingleCPUKernel;

    bpy::class_<prog_factory_t>("ProgramFactory", bpy::no_init)
            .def("create_add_particles", rp::adapt_unique(f_add_particle))
            .def("create_euler_integrator", rp::adapt_unique(f_euler_integrator))
            .def("create_update_forces", rp::adapt_unique(f_forces))
            .def("create_update_neighbor_list", rp::adapt_unique(f_neighbor_list))
            .def("create_reactions_uncontrolled_approximation", rp::adapt_unique(f_reactions_uncontrolled_approximation));

    bpy::class_<add_particle_t>("AddParticle", bpy::init<scpu_kernel_t*>())
            .def("execute", &add_particle_t::execute)
            .def("add_particle", &add_particle_t::addParticle)
            .def("set_particles", &add_particle_t::setParticles);

    bpy::class_<euler_integrator_t>("EulerBDIntegrator", bpy::init<scpu_kernel_t*>())
            .def("execute", &euler_integrator_t::execute);

    bpy::class_<forces_t>("CalculateForces", bpy::init<scpu_kernel_t*>())
            .def("execute", &forces_t::execute);
}