/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * << detailed description >>
 *
 * @file ExportKernelContext.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <pybind11/pybind11.h>
#include <readdy/model/KernelContext.h>
#include <readdy/api/KernelConfiguration.h>

namespace py = pybind11;

using Context = readdy::model::KernelContext;
using KernelConfiguration = readdy::conf::Configuration;

void exportKernelContext(py::module &m) {
    using namespace readdy;
    using namespace py::literals;

    py::class_<Context>(m, "KernelContext")
            .def(py::init<>())
            .def_property("kbt", [](const Context &self) { return self.kBT(); },
                          [](Context &self, scalar kbt) { self.kBT() = kbt; })
            .def("box_volume", &Context::boxVolume)
            .def_property("box_size", [](const Context &self) { return self.boxSize(); },
                          [](Context &self, Context::BoxSize boxSize) { self.boxSize() = boxSize; })
            .def_property("pbc", [](const Context &self) { return self.periodicBoundaryConditions(); },
                          [](Context &self, Context::PeriodicBoundaryConditions pbc) {
                              self.periodicBoundaryConditions() = pbc;
                          })
            .def("bounding_box_vertices", &Context::getBoxBoundingVertices)
            .def("calculate_max_cutoff", &Context::calculateMaxCutoff)
            .def_property("record_reactions_with_positions",
                          [](const Context &self) { return self.recordReactionsWithPositions(); },
                          [](Context &self, bool value) { self.recordReactionsWithPositions() = value; })
            .def_property("record_reaction_counts", [](const Context &self) { return self.recordReactionCounts(); },
                          [](Context &self, bool value) { self.recordReactionCounts() = value; })
            .def("set_kernel_configuration", &Context::setKernelConfiguration)
            ;

}