//
// Created by mho on 8/23/21.
//

#include <pybind11/pybind11.h>
#include "readdy/model/geometry.h"
namespace py = pybind11;

void exportGeometries(py::module &m) {
    using namespace readdy;
    py::class_<model::geometry::Box<scalar>>(m, "Box").def(py::init([](Vec3 v0, Vec3 v1) {
            return model::geometry::Box<scalar>{.v0=v0, .v1=v1};
    }));
    py::class_<model::geometry::Sphere<scalar>>(m, "Sphere").def(py::init([](Vec3 center, scalar radius) {
        return model::geometry::Sphere<scalar>{.center=center, .radius=radius};
    }));
    py::class_<model::geometry::Capsule<scalar>>(m, "Capsule").def(py::init([](Vec3 center, Vec3 direction, scalar radius, scalar length) {
        return model::geometry::Capsule<scalar>{.center=center, .direction=direction / direction.norm(), .radius=radius, .length=length};
    }));
}
