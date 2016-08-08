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

namespace bpy = boost::python;
namespace rp = readdy::py;

using ctx_t = readdy::model::KernelContext;
using particle_t = readdy::model::Particle;

using scpu_nl_t = readdy::kernel::singlecpu::model::SingleCPUNeighborList;
using scpu_nl_box_t = readdy::kernel::singlecpu::model::Box;
using scpu_pd_t = readdy::kernel::singlecpu::model::SingleCPUParticleData;

void exportModelClasses() {
    bpy::class_<std::vector<unsigned long>>("Vec_ulong").def(bpy::vector_indexing_suite<std::vector<unsigned long>>());
    bpy::class_<std::vector<scpu_nl_box_t>>("Vec_box").def(bpy::vector_indexing_suite<std::vector<scpu_nl_box_t>>());

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
            .def("remove_particle", +[](scpu_pd_t &self, particle_t particle) {self.removeParticle(particle);})
            .def("remove_particle", +[](scpu_pd_t &self, std::size_t index) {self.removeParticle(index);})
            .def("is_marked_for_deactivation", &scpu_pd_t::isMarkedForDeactivation)
            .def("get_deactivated_index", &scpu_pd_t::getDeactivatedIndex)
            .def("get_n_deactivated", &scpu_pd_t::getNDeactivated)
            .def("mark_for_deactivation", &scpu_pd_t::markForDeactivation)
            .def("deactivate_marked", &scpu_pd_t::deactivateMarked)
            .def("__getitem__", +[](scpu_pd_t &self, const unsigned int i) {return self[i];});
}