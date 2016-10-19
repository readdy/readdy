//
// Created by mho on 10/08/16.
//

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl_bind.h>

#include <readdy/model/Vec3.h>

#include <readdy/common/logging.h>

namespace bpy = pybind11;

/**
 * Notice: Exporting classes here that are to be shared between prototyping and api module require the base
 * class to use be exported (preferably by the READDY_EXPORT macro defined in common/macros.h).
 */

// module
PYBIND11_PLUGIN (common) {

    if(!readdy::log::console()) {
        spdlog::set_sync_mode();
        auto console = spdlog::stdout_color_mt("console");
        console->set_level(spdlog::level::debug);
        console->set_pattern("[          ] [%Y-%m-%d %H:%M:%S] [%t] [%l] %v");
    }

    bpy::module common("common", "ReaDDy common python module");

    bpy::class_<readdy::model::Vec3>(common, "Vec")
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
            .def("__repr__", [](const readdy::model::Vec3 &self) {
                std::ostringstream stream;
                stream << self;
                return stream.str();
            })
            .def("__getitem__", [](const readdy::model::Vec3 &self, unsigned int i) {
                return self[i];
            });

    return common.ptr();
}