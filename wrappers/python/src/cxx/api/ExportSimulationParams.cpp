/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * @file ExportSimulationParams.cpp
 * @brief << brief description >>
 * @author chrisfroe
 * @date 20.02.19
 */

#include <pybind11/pybind11.h>
//#include <pybind11/stl_bind.h>
//#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <readdy/model/SimulationParams.h>

namespace py = pybind11;

using rvp = py::return_value_policy;
using SimulationParams = readdy::model::SimulationParams;

void exportSimulationParams(py::module &module) {
    //using namespace readdy;
    //using namespace py::literals;

    py::class_<SimulationParams>(module, "SimulationParams")
            .def_readwrite("recordReactionsWithPositions", &SimulationParams::recordReactionsWithPositions)
            .def_readwrite("recordReactionCounts", &SimulationParams::recordReactionCounts)
            .def_readwrite("recordVirial", &SimulationParams::recordVirial)
            .def_readwrite("neighborListCellWidth", &SimulationParams::neighborListCellWidth)
            .def_readwrite("neighborListSkinSize", &SimulationParams::neighborListSkinSize);
}