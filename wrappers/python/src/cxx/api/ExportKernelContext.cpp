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

using KernelContext = readdy::model::KernelContext;
using KernelConfiguration = readdy::conf::Configuration;

void exportKernelContext(py::module &m) {
    using namespace readdy;
    using namespace py::literals;

    py::class_<KernelContext>(m, "KernelContext")
            .def(py::init<>())
            .def_property("kbt", [](const KernelContext &self) { return self.kBT(); },
                          [](KernelContext &self, scalar kbt) { self.kBT() = kbt; })
            .def("box_volume", &KernelContext::boxVolume)
            .def_property("box_size", [](const KernelContext &self) { return self.boxSize(); },
                          [](KernelContext &self, KernelContext::BoxSize boxSize) { self.boxSize() = boxSize; })
            .def_property("pbc", [](const KernelContext &self) { return self.periodicBoundaryConditions(); },
                          [](KernelContext &self, KernelContext::PeriodicBoundaryConditions pbc) {
                              self.periodicBoundaryConditions() = pbc;
                          })
            .def("bounding_box_vertices", &KernelContext::getBoxBoundingVertices)
            .def("calculate_max_cutoff", &KernelContext::calculateMaxCutoff)
            .def_property("record_reactions_with_positions",
                          [](const KernelContext &self) { return self.recordReactionsWithPositions(); },
                          [](KernelContext &self, bool value) { self.recordReactionsWithPositions() = value; })
            .def_property("record_reaction_counts", [](const KernelContext &self) { return self.recordReactionCounts(); },
                          [](KernelContext &self, bool value) { self.recordReactionCounts() = value; })
            .def("set_kernel_configuration", &KernelContext::setKernelConfiguration)
            ;

}