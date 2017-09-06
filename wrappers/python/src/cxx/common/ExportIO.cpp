/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * @file ExportIO.cpp
 * @brief File containing binding points to the readdy.io module
 * @author clonker
 * @date 13.12.16
 * @copyright GNU Lesser General Public License v3.0
 */


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <spdlog/fmt/ostr.h>
#include <h5rd/h5rd.h>

namespace py = pybind11;
using rvp = py::return_value_policy;
namespace io = h5rd;

template<typename T>
void exportDataSet(py::module &io, const std::string &name) {
    using group_t = io::Group;
    using dataset_t = io::DataSet  ;
    std::string base_name = "DataSet_";
    py::class_<dataset_t>(io, (base_name + name).c_str())
            .def("append", [](dataset_t &self, const py::array_t <T> &arr) {
                self.append(io::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            });
}

void exportIO(py::module &io) {
    using namespace pybind11::literals;
    using file_t = io::File;
    using file_node_t = io::Node<file_t>;
    using group_t = io::Group;
    using group_node_t = io::Node<group_t>;

    io.def("unlimited_dims", [] { return io::UNLIMITED_DIMS; });

    py::enum_<file_t::Action>(io, "FileAction")
            .value("CREATE", file_t::Action::CREATE)
            .value("OPEN", file_t::Action::OPEN)
            .export_values();
    py::enum_<file_t::Flag>(io, "FileFlag")
            .value("DEFAULT", file_t::Flag::DEFAULT)
            .value("READ_ONLY", file_t::Flag::READ_ONLY)
            .value("READ_WRITE", file_t::Flag::READ_WRITE)
            .value("OVERWRITE", file_t::Flag::OVERWRITE)
            .value("FAIL_IF_EXISTS", file_t::Flag::FAIL_IF_EXISTS)
            .value("CREATE_NON_EXISTING", file_t::Flag::CREATE_NON_EXISTING)
            .export_values();

    py::class_<file_node_t, std::shared_ptr<file_node_t>>(io, "FileNode")
            .def("create_group", &file_node_t::createGroup, rvp::move)
            .def("subgroups", &file_node_t::subgroups)
            .def("data_sets", &file_node_t::containedDataSets)
            .def("get_subgroup", &file_node_t::getSubgroup);

    py::class_<group_node_t, std::shared_ptr<group_node_t>>(io, "GroupNode")
            .def("create_group", &group_node_t::createGroup, rvp::move)
            .def("subgroups", &group_node_t::subgroups)
            .def("data_sets", &group_node_t::containedDataSets)
            .def("get_subgroup", &group_node_t::getSubgroup);

    py::class_<io::Object, std::shared_ptr<io::Object>>(io, "Object").def("hid", &io::Object::id);

    py::class_<file_t, std::shared_ptr<file_t>, io::Object, file_node_t>(io, "File", py::multiple_inheritance())
            .def_static("open", [](const std::string &path, const io::File::Flag &flag) {
                return io::File::open(path, flag);
            }, "path"_a, "flag"_a = io::File::Flag::READ_WRITE)
            .def_static("create", [](const std::string &path, const io::File::Flag &flag) {
                return io::File::create(path, flag);
            }, "path"_a, "flag"_a = io::File::Flag::OVERWRITE)
            .def("write_short", [](file_t &self, const std::string &name, const py::array_t<short> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_int", [](file_t &self, const std::string &name, const py::array_t<int> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_long", [](file_t &self, const std::string &name, const py::array_t<long> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_float", [](file_t &self, const std::string &name, const py::array_t<float> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_double", [](file_t &self, const std::string &name, const py::array_t<double> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_string", [](file_t &self, const std::string &name, const std::string &data) {
                self.write(name, data);
            })
            .def("flush", &file_t::flush)
            .def("close", &file_t::close)
            .def("create_group", &file_t::createGroup);

    py::class_<group_t, std::shared_ptr<group_t>, io::Object, group_node_t>(io, "Group", py::multiple_inheritance())
            .def("write_short", [](group_t &self, const std::string &name, const py::array_t<short> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_int", [](group_t &self, const std::string &name, const py::array_t<int> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_long", [](group_t &self, const std::string &name, const py::array_t<long> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_float", [](group_t &self, const std::string &name, const py::array_t<float> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_double", [](group_t &self, const std::string &name, const py::array_t<double> &arr) {
                self.write(name, h5rd::dimensions(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_string", [](group_t &self, const std::string &name, const std::string &data) {
                self.write(name, data);
            });

    //exportDataSet<short>(io, "short"); /* DataSet_short */
    //exportDataSet<int>(io, "int"); /* DataSet_int */
    //exportDataSet<long>(io, "long"); /* DataSet_long */
    //exportDataSet<float>(io, "float"); /* DataSet_float */
    //exportDataSet<readdy::scalar>(io, "readdy::scalar"); /* readdy::scalar */

}