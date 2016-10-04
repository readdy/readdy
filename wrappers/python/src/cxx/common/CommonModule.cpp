//
// Created by mho on 10/08/16.
//
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <boost/python.hpp>

#include <vector>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>
#include <boost/uuid/uuid_io.hpp>
#include "../PyConverters.h"

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

namespace bpy = boost::python;

using _rdy_scpu_nl_box_t = readdy::kernel::singlecpu::model::Box;
using vec = readdy::model::Vec3;
using uuid = boost::uuids::uuid;

double vecBracketOperator(vec &self, unsigned int i) {
    return self[i];
}

// module
BOOST_PYTHON_MODULE (common) {
    init_numpy();
    PyEval_InitThreads();

    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    bpy::docstring_options doc_options;
    doc_options.enable_all();

    bpy::class_<vec>("Vec", bpy::init<double, double, double>())
            .def(bpy::self + bpy::self)
            .def(bpy::self - bpy::self)
            .def(double() * bpy::self)
            .def(bpy::self / double())
            .def(bpy::self += bpy::self)
            .def(bpy::self *= double())
            .def(bpy::self == bpy::self)
            .def(bpy::self != bpy::self)
            .def(bpy::self * bpy::self)
            .def(bpy::self_ns::str(bpy::self))
            .def("__getitem__", &vecBracketOperator);

    readdy::py::std_vector_to_python_converter<double>();
    readdy::py::std_pair_to_python_converter<std::vector<double>, std::vector<double>>();
    bpy::class_<std::vector<unsigned long>>("Vec_ulong").def(bpy::vector_indexing_suite<std::vector<unsigned long>>());
    bpy::class_<std::vector<_rdy_scpu_nl_box_t>>("Vec_box").def(
            bpy::vector_indexing_suite<std::vector<_rdy_scpu_nl_box_t>>());
    bpy::class_<std::vector<vec>>("Vecvec").def(boost::python::vector_indexing_suite<std::vector<vec>>());
    bpy::class_<uuid>("uuid", bpy::no_init).def("__str__",
                                                +[](const uuid &uuid) { return boost::uuids::to_string(uuid); });
}