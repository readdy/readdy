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
 * @file Utils.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.07.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <fstream>

#include <spdlog/fmt/ostr.h>

#include <readdy/model/observables/io/TrajectoryEntry.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/IOUtils.h>
#include <readdy/io/BloscFilter.h>

namespace py = pybind11;
using rvp = py::return_value_policy;

using radiusmap = std::map<std::string, readdy::scalar>;

py::tuple convert_readdy_viewer(const std::string &h5name, const std::string &trajName) {
    readdy::log::debug(R"(converting "{}" to readdy viewer format)", h5name);

    readdy::io::BloscFilter bloscFilter;
    bloscFilter.registerFilter();

    auto f = h5rd::File::open(h5name, h5rd::File::Flag::READ_ONLY);

    auto particleInfoH5Type = readdy::model::ioutils::getReactionInfoMemoryType(f->ref());

    // get particle types from config
    std::vector<readdy::model::ioutils::ParticleTypeInfo> types;
    {
        auto config = f->getSubgroup("readdy/config");
        config.read("particle_types", types, &std::get<0>(particleInfoH5Type), &std::get<1>(particleInfoH5Type));
    }

    auto traj = f->getSubgroup("readdy/trajectory/" + trajName);

    // limits
    std::vector<std::size_t> limits;
    traj.read("limits", limits);

    // records
    std::vector<readdy::model::observables::TrajectoryEntry> entries;
    auto trajectoryEntryTypes = readdy::model::observables::util::getTrajectoryEntryTypes(f->ref());
    traj.read("records", entries, &std::get<0>(trajectoryEntryTypes), &std::get<1>(trajectoryEntryTypes));

    auto n_frames = limits.size()/2;
    readdy::log::debug("got n frames: {}", n_frames);

    // map from type name to max number of particles in traj
    std::size_t max_n_particles_per_frame = 0;
    std::vector<std::size_t> shape_nppf {n_frames};
    py::array_t<std::size_t, py::array::c_style> n_particles_per_frame (shape_nppf);
    {
        auto ptr = n_particles_per_frame.mutable_data(0);

        for (std::size_t i = 0; i < limits.size(); i += 2) {
            auto len = limits[i + 1] - limits[i];
            ptr[i/2] = len;
            max_n_particles_per_frame = std::max(max_n_particles_per_frame, len);
        }
    }
    readdy::log::debug("got max n particles: {}", max_n_particles_per_frame);

    py::array_t<double, py::array::c_style> positions_arr {std::vector<std::size_t>{n_frames, max_n_particles_per_frame, 3}};
    py::array_t<std::size_t, py::array::c_style> types_arr {std::vector<std::size_t>{n_frames, max_n_particles_per_frame}};
    py::array_t<std::size_t, py::array::c_style> ids_arr {std::vector<std::size_t>{n_frames, max_n_particles_per_frame}};

    auto positions_ptr = positions_arr.mutable_unchecked();
    auto types_ptr = types_arr.mutable_unchecked();
    auto ids_ptr = ids_arr.mutable_unchecked();

    for (std::size_t i = 0; i < limits.size(); i += 2) {
        auto frame = i/2;
        auto begin = limits[i];
        auto end = limits[i+1];

        std::size_t p = 0;
        for (auto it = entries.begin() + begin; it != entries.begin() + end; ++it, ++p) {
            positions_ptr(frame, p, 0) = it->pos.x;
            positions_ptr(frame, p, 1) = it->pos.y;
            positions_ptr(frame, p, 2) = it->pos.z;
            types_ptr(frame, p) = it->typeId;
            ids_ptr(frame, p) = it->id;
        }
    }

    return py::make_tuple(n_particles_per_frame, positions_arr, types_arr, ids_arr);
}

void
convert_xyz(const std::string &h5name, const std::string &trajName, const std::string &out, bool generateTcl = true, bool tclRuler = false,
            const radiusmap &radii = {}) {
    readdy::log::debug(R"(converting "{}" to "{}")", h5name, out);

    readdy::io::BloscFilter bloscFilter;
    bloscFilter.registerFilter();

    auto f = h5rd::File::open(h5name, h5rd::File::Flag::READ_ONLY);

    auto particleInfoH5Type = readdy::model::ioutils::getReactionInfoMemoryType(f->ref());

    // get particle types from config
    std::vector<readdy::model::ioutils::ParticleTypeInfo> types;
    {
        auto config = f->getSubgroup("readdy/config");
        config.read("particle_types", types, &std::get<0>(particleInfoH5Type), &std::get<1>(particleInfoH5Type));
    }

    // map from type name to max number of particles in traj
    std::unordered_map<readdy::particle_type_type, std::size_t> maxCounts;
    for (const auto &type : types) {
        maxCounts[type.type_id] = 0;
        readdy::log::debug("got type {} with id {} and D {}", type.name, type.type_id, type.diffusion_constant);
    }

    auto traj = f->getSubgroup("readdy/trajectory/" + trajName);

    // limits
    std::vector<std::size_t> limits;
    traj.read("limits", limits);
    // records
    std::vector<readdy::model::observables::TrajectoryEntry> entries;


    auto trajectoryEntryTypes = readdy::model::observables::util::getTrajectoryEntryTypes(f->ref());
    traj.read("records", entries, &std::get<0>(trajectoryEntryTypes), &std::get<1>(trajectoryEntryTypes));

    std::unordered_map<readdy::particle_type_type, std::size_t> typeMapping(types.size());
    {
        std::size_t i = 0;
        for (const auto &type : types) {
            typeMapping[type.type_id] = i++;
        }
    }

    {
        std::vector<std::size_t> currentCounts(types.size());
        for (std::size_t i = 0; i < limits.size(); i += 2) {
            std::fill(currentCounts.begin(), currentCounts.end(), 0);

            auto begin = limits[i];
            auto end = limits[i + 1];

            for (auto it = entries.begin() + begin; it != entries.begin() + end; ++it) {
                currentCounts[typeMapping.at(it->typeId)]++;
            }

            for (const auto &e : typeMapping) {
                maxCounts[e.first] = std::max(currentCounts.at(e.second), maxCounts[e.first]);
            }
        }
    }

    readdy::log::debug("got particle counts:");
    for (const auto &e : maxCounts) {
        readdy::log::debug("\t type id {}: max {} particles", e.first, e.second);
    }

    std::size_t maxParticlesSum = 0;
    {
        for (const auto &e :maxCounts) {
            maxParticlesSum += e.second;
        }
    }

    readdy::log::debug("writing to xyz (n timesteps {})", limits.size() / 2);

    {
        std::fstream fs;
        fs.exceptions(std::fstream::failbit);
        fs.open(out, std::fstream::out | std::fstream::trunc);

        std::vector<std::size_t> currentCounts(types.size());
        std::vector<std::string> xyzPerType(types.size());
        for (std::size_t frame = 0; frame < limits.size(); frame += 2) {

            // number of atoms + comment line (empty)
            fs << maxParticlesSum << std::endl << std::endl;

            std::fill(currentCounts.begin(), currentCounts.end(), 0);
            std::fill(xyzPerType.begin(), xyzPerType.end(), "");

            auto begin = limits[frame];
            auto end = limits[frame + 1];

            for (auto it = entries.begin() + begin; it != entries.begin() + end; ++it) {
                currentCounts[typeMapping.at(it->typeId)]++;
                auto &currentXYZ = xyzPerType.at(typeMapping.at(it->typeId));
                currentXYZ += "type_" + std::to_string(it->typeId) + "\t" + std::to_string(it->pos.x) + "\t" +
                              std::to_string(it->pos.y) + "\t" + std::to_string(it->pos.z) + "\n";
            }

            for (const auto &mappingEntry : typeMapping) {
                auto nGhosts = maxCounts.at(mappingEntry.first) - currentCounts.at(mappingEntry.second);

                fs << xyzPerType.at(mappingEntry.second);
                for (int x = 0; x < nGhosts; ++x) {
                    fs << "type_" + std::to_string(mappingEntry.first) + "\t0\t0\t0\n";
                }
            }
        }

    }

    if (generateTcl) {
        readdy::log::debug("generating tcl script file {}.tcl", out);

        std::fstream fs;
        fs.exceptions(std::fstream::failbit);
        fs.open(out + ".tcl", std::fstream::out | std::fstream::trunc);

        fs << "mol delete top" << std::endl;
        fs << "mol load xyz " << out << std::endl;
        fs << "mol delrep 0 top" << std::endl;
        fs << "display resetview" << std::endl;

        if (tclRuler) {
            fs << "display projection orthographic" << std::endl;
            fs << "ruler grid" << std::endl;
        }

        std::size_t i = 0;
        for (const auto &t : types) {
            auto radius = static_cast<readdy::scalar>(1.0);
            {
                auto it = radii.find(t.name);
                if (it == radii.end()) {
                    readdy::log::warn("type {} had explicitly no specified radius, using 1.0", t.name);
                } else {
                    radius = it ->second;
                }
            }
            fs << "mol representation VDW " << std::to_string(radius * .7) << " 16.0" << std::endl;
            fs << "mol selection name type_" + std::to_string(t.type_id) << std::endl;
            fs << "mol color ColorID " << (i++) << std::endl;
            fs << "mol addrep top" << std::endl;
        }
        fs << "animate goto 0" << std::endl;
        fs << "color Display Background white" << std::endl;
        fs << "molinfo top set {center_matrix} {{{1 0 0 0}{0 1 0 0}{0 0 1 0}{0 0 0 1}}}" << std::endl;
    }

    readdy::log::debug("converting finished");
}

void exportUtils(py::module &m) {
    using namespace pybind11::literals;
    m.def("convert_xyz", &convert_xyz, "h5_file_name"_a, "traj_data_set_name"_a, "xyz_out_file_name"_a,
          "generate_tcl"_a = true, "tcl_with_grid"_a = false, "radii"_a = radiusmap{});
    m.def("convert_readdyviewer", &convert_readdy_viewer, "h5_file_name"_a, "traj_data_set_name"_a);
}