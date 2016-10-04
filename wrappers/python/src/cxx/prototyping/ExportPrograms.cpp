/**
 * << detailed description >>
 *
 * @file Programs.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */
#include <boost/python.hpp>
#include <../PyConverters.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUUpdateNeighborList.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>
#include <../PyFunction.h>
#include "ProgramWrap.h"

namespace bpy = boost::python;
namespace rpy = readdy::py;

using _rdy_particle_t = readdy::model::Particle;
using py_fun_11_t = rpy::PyFunction<_rdy_particle_t(const _rdy_particle_t &)>;
using py_fun_12_t = rpy::PyFunction<void(const _rdy_particle_t &, _rdy_particle_t &, _rdy_particle_t &)>;
using py_fun_21_t = rpy::PyFunction<_rdy_particle_t(const _rdy_particle_t &, const _rdy_particle_t &)>;
using py_fun_22_t = rpy::PyFunction<void(const _rdy_particle_t &, const _rdy_particle_t &, _rdy_particle_t &,
                                         _rdy_particle_t &)>;

using prog_factory_t = readdy::model::programs::ProgramFactory;

using program_t = readdy::model::programs::Program;
using program_wrap_t = readdy::py::Program;

using add_particle_t = readdy::kernel::singlecpu::programs::SingleCPUAddParticleProgram;
using euler_integrator_t = readdy::kernel::singlecpu::programs::SingleCPUEulerBDIntegrator;
using forces_t = readdy::kernel::singlecpu::programs::SingleCPUCalculateForces;
using neighbor_list_t = readdy::kernel::singlecpu::programs::SingleCPUUpdateNeighborList;

using reactions_u_a_t = readdy::kernel::singlecpu::programs::reactions::UncontrolledApproximation;

void exportPrograms() {
    auto f_add_particle = &prog_factory_t::createProgram<add_particle_t>;
    auto f_euler_integrator = &prog_factory_t::createProgram<euler_integrator_t>;
    auto f_forces = &prog_factory_t::createProgram<forces_t>;
    auto f_neighbor_list = &prog_factory_t::createProgram<neighbor_list_t>;
    auto f_reactions_uncontrolled_approximation = &prog_factory_t::createProgram<reactions_u_a_t>;

    using scpu_kernel_t = readdy::kernel::singlecpu::SingleCPUKernel;

    bpy::class_<prog_factory_t>("ProgramFactory", bpy::no_init)
            .def("create_add_particles", rpy::adapt_unique(f_add_particle))
            .def("create_euler_integrator", rpy::adapt_unique(f_euler_integrator))
            .def("create_update_forces", rpy::adapt_unique(f_forces))
            .def("create_update_neighbor_list", rpy::adapt_unique(f_neighbor_list))
            .def("create_reactions_uncontrolled_approximation",
                 rpy::adapt_unique(f_reactions_uncontrolled_approximation));

    bpy::class_<program_wrap_t, boost::noncopyable>("Program").def("execute", bpy::pure_virtual(&program_t::execute));

    bpy::class_<add_particle_t, boost::python::bases<program_t>>("AddParticle", bpy::init<scpu_kernel_t *>())
            .def("execute", &add_particle_t::execute)
            .def("add_particle", &add_particle_t::addParticle)
            .def("set_particles", &add_particle_t::setParticles);

    bpy::class_<euler_integrator_t, boost::python::bases<program_t>>("EulerBDIntegrator", bpy::init<scpu_kernel_t *>())
            .def("execute", &euler_integrator_t::execute);

    bpy::class_<forces_t, boost::python::bases<program_t>>("CalculateForces", bpy::init<scpu_kernel_t *>())
            .def("execute", &forces_t::execute);

    bpy::class_<neighbor_list_t, boost::python::bases<program_t>>("UpdateNeighborList", bpy::init<scpu_kernel_t *>())
            .def("execute", &neighbor_list_t::execute);

    bpy::class_<reactions_u_a_t, boost::python::bases<program_t>>("ReactionsUncontrolledApproximation",
                                                                  bpy::init<scpu_kernel_t *>())
            .def("execute", &reactions_u_a_t::execute)
            .def("register_reaction_scheme_11",
                 +[](reactions_u_a_t &self, const std::string name, bpy::object handler) {
                     self.registerReactionScheme_11(name, py_fun_11_t(handler));
                 })
            .def("register_reaction_scheme_12",
                 +[](reactions_u_a_t &self, const std::string name, bpy::object handler) {
                     self.registerReactionScheme_12(name, py_fun_12_t(handler));
                 })
            .def("register_reaction_scheme_21",
                 +[](reactions_u_a_t &self, const std::string name, bpy::object handler) {
                     self.registerReactionScheme_21(name, py_fun_21_t(handler));
                 })
            .def("register_reaction_scheme_22",
                 +[](reactions_u_a_t &self, const std::string name, bpy::object handler) {
                     self.registerReactionScheme_22(name, py_fun_22_t(handler));
                 });
}