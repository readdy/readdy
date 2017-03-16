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
 * << detailed description >>
 *
 * @file ExportObservables.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#ifndef READDY_MAIN_EXPORTOBSERVABLES_H
#define READDY_MAIN_EXPORTOBSERVABLES_H

#include <pybind11/pybind11.h>
#include <readdy/model/observables/Observables.h>
#include <readdy/api/Simulation.h>
#include <pybind11/numpy.h>
#include <readdy/model/observables/io/Trajectory.h>
#include "PyFunction.h"

namespace py = pybind11;
using rvp = py::return_value_policy;
using sim = readdy::Simulation;
using obs_handle_t = readdy::ObservableHandle;

inline obs_handle_t registerObservable_Reactions(sim &self, unsigned int stride, py::object callback) {
    if (callback.is_none()) {
        return self.registerObservable<readdy::model::observables::Reactions>(stride);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Reactions::result_t)>(callback);
        return self.registerObservable<readdy::model::observables::Reactions>(std::move(pyFun), stride);
    }
}

inline obs_handle_t registerObservable_ReactionCounts(sim &self, unsigned int stride, py::object callback) {
    if (callback.is_none()) {
        return self.registerObservable<readdy::model::observables::ReactionCounts>(stride);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::ReactionCounts::result_t)>(callback);
        return self.registerObservable<readdy::model::observables::ReactionCounts>(std::move(pyFun), stride);
    }
}

inline obs_handle_t
registerObservable_Positions(sim &self, unsigned int stride, pybind11::object callbackFun,
                             std::vector<std::string> types) {
    if (callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::Positions>(stride, types);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Positions::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::Positions>(std::move(pyFun), stride, types);
    }
}

inline obs_handle_t registerObservable_Particles(sim &self, unsigned int stride, pybind11::object callbackFun) {
    if (callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::Particles>(stride);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Particles::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::Particles>(std::move(pyFun), stride);
    }
}


inline obs_handle_t
registerObservable_RadialDistribution(sim &self, unsigned int stride, pybind11::object callbackFun,
                                      py::array_t<double> &binBorders, std::vector<std::string> typeCountFrom,
                                      std::vector<std::string> typeCountTo, double particleToDensity) {
    const auto info = binBorders.request();
    std::vector<double> binBordersVec{};
    binBordersVec.reserve(info.shape[0]);
    const auto data = static_cast<double *>(info.ptr);
    for (auto i = 0; i < info.shape[0]; ++i) binBordersVec.push_back(data[i]);
    if (callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::RadialDistribution>(
                stride, binBordersVec, typeCountFrom, typeCountTo, particleToDensity
        );
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::RadialDistribution::result_t)>(
                callbackFun);
        return self.registerObservable<readdy::model::observables::RadialDistribution>(
                std::move(pyFun), stride, binBordersVec, typeCountFrom, typeCountTo, particleToDensity
        );
    }
}

inline obs_handle_t
registerObservable_CenterOfMass(sim &self, unsigned int stride, const pybind11::object &callbackFun,
                                std::vector<std::string> types) {
    if (callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::CenterOfMass>(stride, types);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::CenterOfMass::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::CenterOfMass>(
                std::move(pyFun), stride, types
        );
    }
}

inline obs_handle_t
registerObservable_HistogramAlongAxisObservable(sim &self, unsigned int stride, const py::object &callbackFun,
                                                py::array_t<double> binBorders, unsigned int axis,
                                                std::vector<std::string> types) {
    const auto info = binBorders.request();
    const auto sizeBorders = info.shape[0];
    auto binBordersData = static_cast<double *>(info.ptr);
    std::vector<double> binBordersVec{};
    binBordersVec.reserve(sizeBorders);
    for (auto i = 0; i < sizeBorders; ++i) {
        binBordersVec.push_back(binBordersData[i]);
    }
    if (callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::HistogramAlongAxis>(
                stride, binBordersVec, types, axis
        );
    } else {
        auto f = readdy::rpy::PyFunction<void(readdy::model::observables::HistogramAlongAxis::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::HistogramAlongAxis>(
                std::move(f), stride, binBordersVec, types, axis
        );
    }
}

inline obs_handle_t
registerObservable_NParticles(sim &self, unsigned int stride, const py::object &callbackFun,
                              std::vector<std::string> types) {
    if (callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::NParticles>(stride, types);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::NParticles::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::NParticles>(std::move(pyFun), stride, types);
    }
}

inline obs_handle_t registerObservable_ForcesObservable(sim &self, unsigned int stride, py::object callbackFun,
                                                 std::vector<std::string> types) {
    if (callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::Forces>(stride, types);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Forces::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::Forces>(std::move(pyFun), stride, types);
    }
}

inline obs_handle_t registerObservable_Trajectory(sim& self, unsigned int stride) {
    return self.registerObservable<readdy::model::observables::Trajectory>(stride);
}

template <typename type_, typename... options>
void exportObservables(py::module &apiModule, py::class_<type_, options...> &simulation) {
    py::class_<obs_handle_t>(apiModule, "ObservableHandle")
            .def(py::init<>())
            .def("enable_write_to_file", &obs_handle_t::enableWriteToFile)
            .def("flush", &obs_handle_t::flush)
            .def("__repr__", [](const obs_handle_t &self) {
                return "ObservableHandle(id=" + std::to_string(self.getId()) + ")";
            });

    using record_t = readdy::model::reactions::ReactionRecord;
    py::class_<record_t>(apiModule, "ReactionRecord")
            .def_readonly("type", &record_t::type)
            .def_readonly("educts", &record_t::educts)
            .def_readonly("products", &record_t::products)
            .def_readonly("types_from", &record_t::types_from)
            .def_readonly("where", &record_t::where)
            .def_readonly("reaction_index", &record_t::reactionIndex)
            .def("__repr__", [](const record_t& self) {
                std::ostringstream ss;
                ss << self;
                return ss.str();
            });

    simulation.def("register_observable_particle_positions", &registerObservable_Positions)
            .def("register_observable_particles", &registerObservable_Particles)
            .def("register_observable_radial_distribution", &registerObservable_RadialDistribution)
            .def("register_observable_histogram_along_axis", &registerObservable_HistogramAlongAxisObservable)
            .def("register_observable_center_of_mass", &registerObservable_CenterOfMass)
            .def("register_observable_n_particles", &registerObservable_NParticles)
            .def("register_observable_forces", &registerObservable_ForcesObservable)
            .def("register_observable_reactions", &registerObservable_Reactions)
            .def("register_observable_reaction_counts", &registerObservable_ReactionCounts)
            .def("register_observable_trajectory", &registerObservable_Trajectory)
            .def("deregister_observable", [](sim &self, const obs_handle_t &handle) {
                self.deregisterObservable(handle.getId());
            });
}
#endif //READDY_MAIN_EXPORTOBSERVABLES_H
