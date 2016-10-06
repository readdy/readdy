//
// Created by mho on 10/08/16.
//
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl_bind.h>

#include "../PybindOpaqueTypes.h"

#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>
#include <boost/uuid/uuid_io.hpp>

namespace bpy = pybind11;

using rdy_scpu_nl_box_t = readdy::kernel::singlecpu::model::Box;
using vec = readdy::model::Vec3;
using uuid = boost::uuids::uuid;

// module
PYBIND11_PLUGIN (common) {

    bpy::module common("common", "ReaDDy common python module");

    bpy::class_<vec>(common, "Vec")
            .def(bpy::init<double, double, double>())
            .def(bpy::self + bpy::self)
            .def(bpy::self - bpy::self)
            .def(double() * bpy::self)
            .def(bpy::self / double())
            .def(bpy::self += bpy::self)
            .def(bpy::self *= double())
            .def(bpy::self == bpy::self)
            .def(bpy::self != bpy::self)
            .def(bpy::self * bpy::self)
            .def("__repr__", [](const vec &self) {
                std::ostringstream stream;
                stream << self;
                return stream.str();
            })
            .def("__getitem__", [](const vec &self, unsigned int i) {
                return self[i];
            });

    bpy::class_<uuid>(common, "uuid").def("__str__", [](const uuid &uuid) { return boost::uuids::to_string(uuid); });

    bpy::bind_vector<std::vector<double>>(common, "StdVectorDouble");
    bpy::bind_vector<std::vector<unsigned long>>(common, "StdVectorUnsignedLong");
    bpy::bind_vector<std::vector<readdy::kernel::singlecpu::model::Box>>(common, "StdVectorSingleCPUBox");
    bpy::bind_vector<std::vector<readdy::model::Vec3>>(common, "StdVectorVec3");

    return common.ptr();
}