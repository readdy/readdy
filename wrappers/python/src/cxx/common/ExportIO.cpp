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
#include <readdy/io/File.h>
#include <readdy/io/DataSet.h>

namespace py = pybind11;
using rvp = py::return_value_policy;
namespace io = readdy::io;

template<typename T>
void exportDataSet(py::module &io, const std::string &name) {
    using group_t = io::Group;
    using dataset_t = io::DataSet<T, false>;
    std::string base_name = "DataSet_";
    py::class_<dataset_t>(io, (base_name + name).c_str())
            .def(py::init<const std::string &, const group_t &, const std::vector<io::h5::dims_t> &,
                    const std::vector<io::h5::dims_t> &>())
            .def("append", [](dataset_t &self, const py::array_t <T> &arr) {
                self.append(std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            });
}

void exportIO(py::module &io) {
    using file_t = io::File;
    using group_t = io::Group;

    io.def("unlimited_dims", [] { return io::h5::UNLIMITED_DIMS; });

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

    py::class_<file_t>(io, "File")
            .def(py::init<const std::string &, file_t::Action, file_t::Flag>(), py::arg("path"), py::arg("action"),
                 py::arg("flag") = file_t::Flag::OVERWRITE)
            .def(py::init<const std::string &, file_t::Action, const std::vector<file_t::Flag> &>())
            .def("write_short", [](file_t &self, const std::string &name, const py::array_t<short> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_int", [](file_t &self, const std::string &name, const py::array_t<int> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_long", [](file_t &self, const std::string &name, const py::array_t<long> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_float", [](file_t &self, const std::string &name, const py::array_t<float> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_double", [](file_t &self, const std::string &name, const py::array_t<double> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_string", [](file_t &self, const std::string &name, const std::string &data) {
                self.write(name, data);
            })
            .def("flush", &file_t::flush)
            .def("close", &file_t::close)
            .def("create_group", &file_t::createGroup, rvp::move)
            .def("get_root_group", &file_t::getRootGroup, rvp::reference_internal);

    py::class_<group_t>(io, "Group")
            .def("write_short", [](group_t &self, const std::string &name, const py::array_t<short> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_int", [](group_t &self, const std::string &name, const py::array_t<int> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_long", [](group_t &self, const std::string &name, const py::array_t<long> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_float", [](group_t &self, const std::string &name, const py::array_t<float> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_double", [](group_t &self, const std::string &name, const py::array_t<double> &arr) {
                self.write(name, std::vector<io::h5::dims_t>(arr.shape(), arr.shape() + arr.ndim()), arr.data());
            })
            .def("write_string", [](group_t &self, const std::string &name, const std::string &data) {
                self.write(name, data);
            })
            .def("create_group", &group_t::createGroup, rvp::move)
            .def("subgroups", &group_t::subgroups)
            .def("data_sets", &group_t::contained_data_sets)
            .def("get_subgroup", &group_t::subgroup);

    exportDataSet<short>(io, "short"); /* DataSet_short */
    exportDataSet<int>(io, "int"); /* DataSet_int */
    exportDataSet<long>(io, "long"); /* DataSet_long */
    exportDataSet<float>(io, "float"); /* DataSet_float */
    exportDataSet<double>(io, "double"); /* DataSet_double */

}