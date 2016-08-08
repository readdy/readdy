/**
 * << detailed description >>
 *
 * @file Model.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 08.08.16
 */

#include <boost/python.hpp>

#include <readdy/model/Particle.h>
#include <PyConverters.h>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <readdy/kernel/singlecpu/SingleCPUKernelStateModel.h>
#include "ModelWrap.h"

namespace bpy = boost::python;
namespace mpl = boost::mpl;
namespace rp = readdy::py;

using uuid_t = boost::uuids::uuid;

using ctx_t = readdy::model::KernelContext;
using particle_t = readdy::model::Particle;
using vec_t = readdy::model::Vec3;

using scpu_model_t = readdy::kernel::singlecpu::SingleCPUKernelStateModel;
using scpu_model_wrap_t = readdy::py::Model;

using scpu_nl_t = readdy::kernel::singlecpu::model::SingleCPUNeighborList;
using scpu_nl_box_t = readdy::kernel::singlecpu::model::Box;
using scpu_pd_t = readdy::kernel::singlecpu::model::SingleCPUParticleData;

using default_policy = bpy::objects::default_iterator_call_policies;

void exportModelClasses() {

    bpy::class_<std::vector<unsigned long>>("Vec_ulong").def(bpy::vector_indexing_suite<std::vector<unsigned long>>());
    bpy::class_<std::vector<scpu_nl_box_t>>("Vec_box").def(bpy::vector_indexing_suite<std::vector<scpu_nl_box_t>>());

    bpy::class_<particle_t>("Particle", boost::python::init<double, double, double, unsigned int>())
            .add_property("pos", +[](particle_t &self) { return self.getPos(); }, &particle_t::setPos)
            .add_property("type", &particle_t::getType, &particle_t::setType)
            .add_property("id", bpy::make_function(&particle_t::getId, bpy::return_value_policy<bpy::reference_existing_object>()))
            .def(bpy::self == bpy::self)
            .def(bpy::self != bpy::self);

    bpy::class_<ctx_t, boost::noncopyable>("Context", bpy::no_init)
            .add_property("kbt", &ctx_t::getKBT, &ctx_t::setKBT)
            .add_property("box_size",
                          +[](ctx_t &self) { return readdy::model::Vec3(self.getBoxSize()); },
                          +[](ctx_t &self, readdy::model::Vec3 vec) { self.setBoxSize(vec[0], vec[1], vec[2]); })
            .def("set_diffusion_constant", &ctx_t::setDiffusionConstant);

    bpy::class_<scpu_model_wrap_t, boost::noncopyable>("Model", bpy::init<ctx_t *>())
            .def("remove_particle", &scpu_model_t::removeParticle, &scpu_model_wrap_t::default_removeParticle)
            .def("get_particle_positions", &scpu_model_t::getParticlePositions, &scpu_model_wrap_t::default_getParticlePositions)
            .def("get_energy", &scpu_model_t::getEnergy, &scpu_model_wrap_t::default_getEnergy)
            .def("increase_energy", &scpu_model_t::increaseEnergy, &scpu_model_wrap_t::default_increaseEnergy)
            .def("get_particle_data", &scpu_model_t::getParticleData, &scpu_model_wrap_t::default_getParticleData, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("get_neighbor_list", &scpu_model_t::getNeighborList, &scpu_model_wrap_t::default_getNeighborList, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("get_particles", &scpu_model_t::getParticles, &scpu_model_wrap_t::default_getParticles);

    bpy::class_<scpu_nl_t, boost::noncopyable>("NeighborList", bpy::init<ctx_t *>())
            .def("create", &scpu_nl_t::create)
            .def("setup_neighboring_boxes", &scpu_nl_t::setupNeighboringBoxes)
            .def("setup_boxes", &scpu_nl_t::setupBoxes)
            .def("fill_boxes", &scpu_nl_t::fillBoxes);
    bpy::class_<scpu_nl_box_t>("NeighborListBox", bpy::init<long, long, long, long>())
            .def("add_neighbor", &scpu_nl_box_t::addNeighbor)
            .def_readonly("i", &scpu_nl_box_t::i)
            .def_readonly("j", &scpu_nl_box_t::j)
            .def_readonly("k", &scpu_nl_box_t::k)
            .def_readonly("id", &scpu_nl_box_t::id)
            .def_readwrite("particle_indices", &scpu_nl_box_t::particleIndices)
            .def_readwrite("neighboring_boxes", &scpu_nl_box_t::neighboringBoxes);

    bpy::class_<scpu_pd_t, boost::noncopyable>("ParticleData")
            .def("swap", &scpu_pd_t::swap)
            .def("size", &scpu_pd_t::size)
            .def("max_size", &scpu_pd_t::max_size)
            .def("empty", &scpu_pd_t::empty)
            .def("clear", &scpu_pd_t::clear)
            .def("add_particle", &scpu_pd_t::addParticle)
            .def("add_particles", &scpu_pd_t::addParticles)
            .def("remove_particle", +[](scpu_pd_t &self, particle_t particle) { self.removeParticle(particle); })
            .def("remove_particle", +[](scpu_pd_t &self, std::size_t index) { self.removeParticle(index); })
            .def("is_marked_for_deactivation", &scpu_pd_t::isMarkedForDeactivation)
            .def("get_deactivated_index", &scpu_pd_t::getDeactivatedIndex)
            .def("get_n_deactivated", &scpu_pd_t::getNDeactivated)
            .def("mark_for_deactivation", &scpu_pd_t::markForDeactivation)
            .def("deactivate_marked", &scpu_pd_t::deactivateMarked)
            .add_property("ids", bpy::range(&scpu_pd_t::cbegin_ids, &scpu_pd_t::cend_ids))
            .add_property("positions", bpy::range(&scpu_pd_t::cbegin_positions, &scpu_pd_t::cend_positions))
            .add_property("forces", bpy::range(&scpu_pd_t::cbegin_forces, &scpu_pd_t::cend_forces))
            .add_property("types", bpy::range(&scpu_pd_t::cbegin_types, &scpu_pd_t::cend_types))
            .def("__getitem__", +[](scpu_pd_t &self, const unsigned int i) { return self[i]; });
}