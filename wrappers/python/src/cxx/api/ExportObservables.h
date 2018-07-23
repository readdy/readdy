/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * << detailed description >>
 *
 * @file ExportObservables.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GPL-3
 */

#ifndef READDY_MAIN_EXPORTOBSERVABLES_H
#define READDY_MAIN_EXPORTOBSERVABLES_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <readdy/model/observables/Observables.h>
#include <readdy/api/Simulation.h>
#include <readdy/model/observables/io/Trajectory.h>
#include "PyFunction.h"
#include "../common/ReadableReactionRecord.h"

namespace py = pybind11;
using rvp = py::return_value_policy;
using sim = readdy::Simulation;
using obs_handle_t = readdy::ObservableHandle;

inline obs_handle_t registerObservable_Reactions(sim &self, unsigned int stride,
                                                 const py::object& callback = py::none()) {
    self.context().recordReactionsWithPositions() = true;
    auto obs = self.observe().reactions(stride);
    if (callback.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto internalCallback = [&self, callback](const readdy::model::observables::Reactions::result_type &reactions) mutable {
            py::gil_scoped_acquire gil;
            std::vector<rpy::ReadableReactionRecord> converted {};
            converted.reserve(reactions.size());
            const auto &reactionRegistry = self.context().reactions();
            for(const auto &reaction : reactions) {
                auto name = reactionRegistry.nameOf(reaction.id);
                converted.emplace_back(rpy::convert(reaction, name));
            }
            for(const auto &c : converted) {
                std::stringstream ss;
                ss << c;

            }
            callback(converted);
        };
        return self.registerObservable(std::move(obs), internalCallback);
    }
}

inline obs_handle_t registerObservable_Topologies(sim &self, readdy::stride_type stride,
                                                  const py::object &callback = py::none()) {
    auto obs = self.observe().topologies(stride);
    if(callback.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Topologies::result_type)>(callback);
        return self.registerObservable(std::move(obs), pyFun);
    }
}

inline obs_handle_t registerObservable_ReactionCounts(sim &self, unsigned int stride,
                                                      const py::object& callback = py::none()) {
    self.context().recordReactionCounts() = true;
    auto obs = self.observe().reactionCounts(stride);
    if (callback.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto internalCallback = [&self, callback](const readdy::model::observables::ReactionCounts::result_type &counts) mutable {
            py::gil_scoped_acquire gil;
            std::unordered_map<std::string, std::size_t> converted;
            const auto &reactionRegistry = self.context().reactions();
            for(const auto &e : counts) {
                converted[reactionRegistry.nameOf(e.first)] = e.second;
            }
            callback(converted);
        };
        return self.registerObservable(std::move(obs), internalCallback);
    }
}

inline obs_handle_t
registerObservable_Positions(sim &self, unsigned int stride,
                             const std::vector<std::string> &types, const pybind11::object& callbackFun = py::none()) {
    auto obs = self.observe().positions(stride, types);
    if (callbackFun.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Positions::result_type)>(callbackFun);
        return self.registerObservable(std::move(obs), pyFun);
    }
}

inline obs_handle_t registerObservable_Particles(sim &self, unsigned int stride,
                                                 const pybind11::object& callbackFun = py::none()) {
    auto obs = self.observe().particles(stride);
    if (callbackFun.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto internalCallback = [&self, callbackFun](const readdy::model::observables::Particles::result_type &r) mutable {
            using particle_type = std::string;
            using particle_id_type = readdy::model::Particle::id_type;
            using result_type = std::tuple<std::vector<particle_type>, std::vector<particle_id_type>, std::vector<readdy::Vec3>>;
            py::gil_scoped_acquire gil;
            result_type result;
            std::get<0>(result).reserve(std::get<0>(r).size());
            std::get<1>(result) = std::get<1>(r);
            std::get<2>(result) = std::get<2>(r);

            auto &names = std::get<0>(result);
            const auto &types = self.context().particleTypes();

            for(const auto particleType : std::get<0>(r)) {
                names.push_back(types.nameOf(particleType));
            }
            callbackFun(result);
        };
        return self.registerObservable(std::move(obs), internalCallback);
    }
}


inline obs_handle_t
registerObservable_RadialDistribution(sim &self, unsigned int stride, py::array_t<readdy::scalar> &binBorders,
                                      const std::vector<std::string> &typeCountFrom, const std::vector<std::string> &typeCountTo,
                                      readdy::scalar particleToDensity, const pybind11::object& callbackFun = py::none()) {
    const auto info = binBorders.request();
    std::vector<readdy::scalar> binBordersVec{};
    binBordersVec.reserve(static_cast<std::size_t>(info.shape[0]));
    const auto data = static_cast<readdy::scalar *>(info.ptr);
    for (auto i = 0; i < info.shape[0]; ++i) binBordersVec.push_back(data[i]);
    auto obs = self.observe().radialDistribution(stride, binBordersVec, typeCountFrom, typeCountTo, particleToDensity);
    if (callbackFun.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::RadialDistribution::result_type)>(
                callbackFun);
        return self.registerObservable(std::move(obs), pyFun);
    }
}

inline obs_handle_t
registerObservable_HistogramAlongAxisObservable(sim &self, unsigned int stride, py::array_t<readdy::scalar> binBorders,
                                                unsigned int axis, std::vector<std::string> types,
                                                const py::object &callbackFun = py::none()) {
    const auto info = binBorders.request();
    const auto sizeBorders = info.shape[0];
    auto binBordersData = static_cast<readdy::scalar *>(info.ptr);
    std::vector<readdy::scalar> binBordersVec{};
    binBordersVec.reserve((unsigned long) sizeBorders);
    for (auto i = 0; i < sizeBorders; ++i) {
        binBordersVec.push_back(binBordersData[i]);
    }
    auto obs = self.observe().histogramAlongAxis(stride, binBordersVec, std::move(types), axis);
    if (callbackFun.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto f = readdy::rpy::PyFunction<void(readdy::model::observables::HistogramAlongAxis::result_type)>(callbackFun);
        return self.registerObservable(std::move(obs), f);
    }
}

inline obs_handle_t
registerObservable_NParticles(sim &self, unsigned int stride, std::vector<std::string> types,
                              const py::object &callbackFun = py::none()) {
    auto obs = self.observe().nParticles(stride, std::move(types));
    if (callbackFun.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::NParticles::result_type)>(callbackFun);
        return self.registerObservable(std::move(obs), pyFun);
    }
}

inline obs_handle_t registerObservable_ForcesObservable(sim &self, unsigned int stride, std::vector<std::string> types,
                                                        const py::object& callbackFun = py::none()) {
    auto obs = self.observe().forces(stride, std::move(types));
    if (callbackFun.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Forces::result_type)>(callbackFun);
        return self.registerObservable(std::move(obs), pyFun);
    }
}

inline obs_handle_t registerObservable_Energy(sim &self, unsigned int stride, const py::object &callbackFun = py::none()) {
    auto obs = self.observe().energy(stride);
    if(callbackFun.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Energy::result_type)>(callbackFun);
        return self.registerObservable(std::move(obs), pyFun);
    }
}

inline obs_handle_t registerObservable_Virial(sim &self, readdy::stride_type stride,
                                              const py::object &callback = py::none()) {
    auto obs = self.observe().virial(stride);
    if(callback.is_none()) {
        return self.registerObservable(std::move(obs));
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Virial::result_type)>(callback);
        return self.registerObservable(std::move(obs), pyFun);
    }
}

inline obs_handle_t registerObservable_Trajectory(sim& self, unsigned int stride) {
    return self.registerObservable(self.observe().trajectory(stride));
}

inline obs_handle_t registerObservable_FlatTrajectory(sim& self, unsigned int stride) {
    return self.registerObservable(self.observe().flatTrajectory(stride));
}

template <typename type_, typename... options>
void exportObservables(py::module &apiModule, py::class_<type_, options...> &simulation) {
    using namespace pybind11::literals;
    py::class_<obs_handle_t>(apiModule, "ObservableHandle")
            .def("enable_write_to_file", &obs_handle_t::enableWriteToFile, "file"_a, "data_set_name"_a, "chunk_size"_a)
            .def("flush", &obs_handle_t::flush)
            .def("__repr__", [](const obs_handle_t &self) {
                return "ObservableHandle(type=" + self.type() + ")";
            });

    using record_t = readdy::model::reactions::ReactionRecord;
    py::class_<record_t>(apiModule, "ReactionRecord")
            .def_property_readonly("type", [](const record_t &self) { return self.type; })
            .def_property_readonly("educts", [](const record_t &self) { return self.educts; })
            .def_property_readonly("products", [](const record_t &self) { return self.products; })
            .def_property_readonly("types_from", [](const record_t &self) { return self.types_from; })
            .def_property_readonly("where", [](const record_t &self) { return self.where; })
            .def_property_readonly("id", [](const record_t &self) { return self.id; })
            .def("__repr__", [](const record_t& self) {
                std::ostringstream ss;
                ss << self;
                return ss.str();
            });

    using TopologyRecord = readdy::model::top::TopologyRecord;
    py::class_<TopologyRecord>(apiModule, "TopologyRecord")
            .def_property_readonly("particles", [](const TopologyRecord &self) {return self.particleIndices;})
            .def_property_readonly("edges", [](const TopologyRecord &self) {return self.edges; });

    simulation.def("register_observable_particle_positions", &registerObservable_Positions,
                   "stride"_a, "types"_a, "callback"_a = py::none())
            .def("register_observable_particles", &registerObservable_Particles, "stride"_a, "callback"_a = py::none())
            .def("register_observable_radial_distribution", &registerObservable_RadialDistribution,
                 "stride"_a, "bin_borders"_a, "type_count_from"_a, "type_count_to"_a,
                 "particle_to_density"_a, "callback"_a = py::none())
            .def("register_observable_histogram_along_axis", &registerObservable_HistogramAlongAxisObservable,
                 "stride"_a, "bin_borders"_a, "axis"_a, "types"_a, "callback"_a = py::none())
            .def("register_observable_n_particles", &registerObservable_NParticles,
                 "stride"_a, "types"_a, "callback"_a = py::none())
            .def("register_observable_forces", &registerObservable_ForcesObservable,
                 "stride"_a, "types"_a, "callback"_a = py::none())
            .def("register_observable_energy", &registerObservable_Energy, "stride"_a, "callback"_a = py::none())
            .def("register_observable_reactions", &registerObservable_Reactions, "stride"_a, "callback"_a = py::none())
            .def("register_observable_reaction_counts", &registerObservable_ReactionCounts,
                 "stride"_a, "callback"_a = py::none())
            .def("register_observable_trajectory", &registerObservable_Trajectory, "stride"_a)
            .def("register_observable_flat_trajectory", &registerObservable_FlatTrajectory, "stride"_a)
            .def("register_observable_virial", &registerObservable_Virial, "stride"_a, "callback"_a=py::none())
            .def("register_observable_topologies", &registerObservable_Topologies, "stride"_a, "callback"_a=py::none());
}
#endif //READDY_MAIN_EXPORTOBSERVABLES_H
