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
 * @file Utils.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.07.17
 * @copyright BSD-3
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>

#include <fstream>
#include <utility>

#include <spdlog/fmt/ostr.h>

#include <readdy/model/observables/io/TrajectoryEntry.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/IOUtils.h>
#include <readdy/io/BloscFilter.h>
#include <readdy/model/reactions/ReactionRecord.h>
#include <readdy/model/topologies/TopologyRecord.h>
#include "ReadableReactionRecord.h"

namespace py = pybind11;
using rvp = py::return_value_policy;

using radiusmap = std::map<std::string, readdy::scalar>;

py::tuple convert_readdy_viewer(const std::string &h5name, const std::string &trajName, std::size_t from,
                                std::size_t to, std::size_t stride) {
    readdy::log::debug(R"(converting "{}" to readdy viewer format)", h5name);

    readdy::io::BloscFilter bloscFilter;
    bloscFilter.registerFilter();

    auto f = h5rd::File::open(h5name, h5rd::File::Flag::READ_ONLY);

    auto particleInfoH5Type = readdy::model::ioutils::getParticleTypeInfoType(f->ref());

    // get particle types from config
    std::vector<readdy::model::ioutils::ParticleTypeInfo> types;
    {
        auto config = f->getSubgroup("readdy/config");
        config.read("particle_types", types, &std::get<0>(particleInfoH5Type), &std::get<1>(particleInfoH5Type));
    }

    auto traj = f->getSubgroup("readdy/trajectory/" + trajName);

    // limits
    std::vector<std::size_t> limits;
    {
        if (stride > 1) {
            traj.read("limits", limits, {stride, 1});
        } else {
            traj.read("limits", limits);
        }
        auto n = limits.size() / 2;

        from = std::min(n, from);
        to = std::min(n, to);

        if (from == to) {
            throw std::invalid_argument(fmt::format("not enough frames to cover range ({}, {}]", from, to));
        } else {
            limits = std::vector<std::size_t>(limits.begin() + 2 * from, limits.begin() + 2 * to);
        }
    }

    auto n_frames = limits.size() / 2;
    readdy::log::debug("got n frames: {}", n_frames);

    if (n_frames > to - from) {
        throw std::logic_error("something major went wrong here");
    }

    // map from type name to max number of particles in traj
    std::size_t max_n_particles_per_frame = 0;
    std::vector<std::size_t> shape_nppf{n_frames};
    py::array_t<std::size_t, py::array::c_style> n_particles_per_frame(shape_nppf);
    {
        auto ptr = n_particles_per_frame.mutable_data(0);

        for (std::size_t i = 0; i < limits.size(); i += 2) {
            auto len = limits[i + 1] - limits[i];
            ptr[i / 2] = len;
            max_n_particles_per_frame = std::max(max_n_particles_per_frame, len);
        }
    }
    readdy::log::debug("got max n particles: {}", max_n_particles_per_frame);

    py::array_t<double, py::array::c_style> positions_arr{
            std::vector<std::size_t>{n_frames, max_n_particles_per_frame, 3}};
    py::array_t<std::size_t, py::array::c_style> types_arr{
            std::vector<std::size_t>{n_frames, max_n_particles_per_frame}};
    py::array_t<std::size_t, py::array::c_style> ids_arr{std::vector<std::size_t>{n_frames, max_n_particles_per_frame}};

    auto positions_ptr = positions_arr.mutable_unchecked();
    auto types_ptr = types_arr.mutable_unchecked();
    auto ids_ptr = ids_arr.mutable_unchecked();

    auto trajectoryEntryTypes = readdy::model::observables::util::getTrajectoryEntryTypes(f->ref());

    for (std::size_t i = 0; i < limits.size(); i += 2) {
        auto frame = i / 2;
        auto begin = limits[i];
        auto end = limits[i + 1];

        // records
        std::vector<readdy::model::observables::TrajectoryEntry> entries;
        traj.readSelection("records", entries, &std::get<0>(trajectoryEntryTypes), &std::get<1>(trajectoryEntryTypes),
                           {begin}, {1}, {end - begin});

        if (entries.size() != end - begin) {
            throw std::runtime_error(fmt::format("this should not happen, the flat selection was {} to {} "
                                                 "but we got {} entries", begin, end, entries.size()));
        }

        std::size_t p = 0;
        for (auto it = entries.begin(); it != entries.end(); ++it, ++p) {
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
convert_xyz(const std::string &h5name, const std::string &trajName, const std::string &out, bool generateTcl = true,
            bool tclRuler = false, const radiusmap &radii = {},
            const std::unordered_map<std::string, unsigned int> &colorIds = {},
            const std::array<readdy::scalar, 3> boxSize = std::array<readdy::scalar, 3>{{0.,0.,0.}}) {
    readdy::log::debug(R"(converting "{}" to "{}")", h5name, out);

    readdy::io::BloscFilter bloscFilter;
    bloscFilter.registerFilter();

    auto f = h5rd::File::open(h5name, h5rd::File::Flag::READ_ONLY);

    auto particleInfoH5Type = readdy::model::ioutils::getParticleTypeInfoType(f->ref());

    // get particle types from config
    std::vector<readdy::model::ioutils::ParticleTypeInfo> types;
    {
        auto config = f->getSubgroup("readdy/config");
        config.read("particle_types", types, &std::get<0>(particleInfoH5Type), &std::get<1>(particleInfoH5Type));
    }

    // map from type name to max number of particles in traj
    std::unordered_map<readdy::ParticleTypeId, std::size_t> maxCounts;
    for (const auto &type : types) {
        maxCounts[type.type_id] = 0;
        readdy::log::debug("got type {} with id {} and D ({}, {}, {})", type.name, type.type_id,
                           type.diffusion_constant[0], type.diffusion_constant[1], type.diffusion_constant[2]);
    }

    auto traj = f->getSubgroup("readdy/trajectory/" + trajName);

    // limits
    std::vector<std::size_t> limits;
    traj.read("limits", limits);
    // records
    std::vector<readdy::model::observables::TrajectoryEntry> entries;


    auto trajectoryEntryTypes = readdy::model::observables::util::getTrajectoryEntryTypes(f->ref());
    traj.read("records", entries, &std::get<0>(trajectoryEntryTypes), &std::get<1>(trajectoryEntryTypes));

    std::unordered_map<readdy::ParticleTypeId, std::size_t> typeMapping(types.size());
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
        std::unordered_set<unsigned int> predefinedColorIds;
        {
            predefinedColorIds.reserve(colorIds.size());
            std::for_each(colorIds.begin(), colorIds.end(), [&predefinedColorIds](const auto &item) {
                predefinedColorIds.emplace(std::get<1>(item));
            });
        }
        for (const auto &t : types) {
            auto radius = static_cast<readdy::scalar>(1.0);
            {
                auto it = radii.find(t.name);
                if (it == radii.end()) {
                    readdy::log::warn("type {} had explicitly no specified radius, using 1.0", t.name);
                } else {
                    radius = it->second;
                }
            }
            std::size_t colorId = i;
            {
                auto it = colorIds.find(t.name);
                if(it != colorIds.end()) {
                    colorId = it->second;
                } else {
                    while(predefinedColorIds.find(colorId) != predefinedColorIds.end()) {
                        ++colorId;
                        ++i;
                    }
                }
            }
            fs << "mol representation VDW " << std::to_string(radius * .7) << " 16.0" << std::endl;
            fs << "mol selection name type_" + std::to_string(t.type_id) << std::endl;
            fs << "mol color ColorID " << colorId << std::endl;
            fs << "mol addrep top" << std::endl;
            i++;
        }
        fs << "animate goto 0" << std::endl;
        fs << "color Display Background white" << std::endl;
        fs << "molinfo top set {center_matrix} {{{1 0 0 0}{0 1 0 0}{0 0 1 0}{0 0 0 1}}}" << std::endl;

        if (not (boxSize.at(0) == 0. and boxSize.at(1) == 0. and boxSize.at(2) == 0.)) {
            readdy::log::debug("drawing simulation box of size [{} {} {}]", boxSize.at(0), boxSize.at(1), boxSize.at(2));

            fs << "set x " << - boxSize.at(0) / 2. << std::endl;
            fs << "set y " << - boxSize.at(1) / 2. << std::endl;
            fs << "set z " << - boxSize.at(2) / 2. << std::endl;
            fs << "set dx " << boxSize.at(0) << std::endl;
            fs << "set dy " << boxSize.at(1) << std::endl;
            fs << "set dz " << boxSize.at(2) << std::endl;
            fs << "set xup [expr $x + $dx]" << std::endl;
            fs << "set yup [expr $y + $dy]" << std::endl;
            fs << "set zup [expr $z + $dz]" << std::endl;
            fs << "set s solid" << std::endl;
            fs << "set w 2" << std::endl;
            fs << "set c 16" << std::endl;
            fs << "graphics top color $c" << std::endl;
            fs << R"(graphics top line "$x $y $z" "$x $y $zup" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$x $y $z" "$x $yup $z" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$x $y $z" "$xup $y $z" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$xup $y $zup" "$x $y $zup" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$xup $y $zup" "$xup $y $z" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$xup $y $zup" "$xup $yup $zup" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$x  $yup  $zup" "$x $y $zup" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$x  $yup  $zup" " $xup  $yup  $zup" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$x  $yup  $zup" "$x  $yup $z" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$xup  $yup $z" "$x  $yup $z" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$xup  $yup $z" " $xup  $yup  $zup" width $w style $s)" << std::endl;
            fs << R"(graphics top line "$xup  $yup $z" " $xup $y $z" width $w style $s)" << std::endl;
        }
    }

    readdy::log::debug("converting finished");
}

std::vector<std::vector<rpy::ReadableReactionRecord>>
read_reactions_obs(const std::string &filename, const std::string &name) {
    readdy::io::BloscFilter bloscFilter;
    bloscFilter.registerFilter();

    auto f = h5rd::File::open(filename, h5rd::File::Flag::READ_ONLY);

    auto reactionInfoH5Type = readdy::model::ioutils::getReactionInfoMemoryType(f->ref());

    // get reaction info from config
    std::unordered_map<readdy::ReactionId, readdy::model::ioutils::ReactionInfo> reactionsMap;
    {
        std::vector<readdy::model::ioutils::ReactionInfo> reactionInfo;
        auto config = f->getSubgroup("readdy/config/");
        config.read("registered_reactions", reactionInfo, &std::get<0>(reactionInfoH5Type),
                    &std::get<1>(reactionInfoH5Type));
        for (const auto &info : reactionInfo) {
            reactionsMap[info.id] = info;
        }
    }

    std::vector<std::vector<readdy::model::reactions::ReactionRecord>> records;
    {
        auto reactionRecordH5Types = readdy::model::observables::util::getReactionRecordTypes(f->ref());
        auto dsgroup = f->getSubgroup("readdy/observables/" + name);
        dsgroup.readVLEN("records", records, &std::get<0>(reactionRecordH5Types), &std::get<1>(reactionRecordH5Types));
    }

    std::vector<std::vector<rpy::ReadableReactionRecord>> result;
    result.reserve(records.size());
    for (const auto &reactions : records) {
        result.emplace_back();
        auto &readableRecords = result.back();
        readableRecords.reserve(reactions.size());

        for (const auto &reaction : reactions) {
            readableRecords.push_back(rpy::convert(reaction, reactionsMap.at(reaction.id).name));
        }

    }

    return result;
}

struct TrajectoryParticle {
    TrajectoryParticle(std::string type, std::string flavor, const std::array<readdy::scalar, 3> &pos,
                       readdy::ParticleId id, readdy::TimeStep t)
            : type(std::move(type)), flavor(std::move(flavor)), position(pos), id(id), t(t) {}

    std::string type;
    std::string flavor;
    std::array<readdy::scalar, 3> position;
    readdy::ParticleId id;
    readdy::TimeStep t;
};

std::string repr(const TrajectoryParticle &p) {
    std::stringstream ss;

    ss << "Particle[id=" << p.id << ", type=" << p.type << ", time=" << p.t << ", flavor=" << p.flavor
       << ", position=(" << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << ")]";

    return ss.str();
}

using TopologyRecord = readdy::model::top::TopologyRecord;

std::tuple<std::vector<readdy::TimeStep>, std::vector<std::vector<TopologyRecord>>>
readTopologies(const std::string &filename, const std::string &groupName, std::size_t from,
               std::size_t to, std::size_t stride) {
    readdy::io::BloscFilter bloscFilter;
    bloscFilter.registerFilter();

    std::vector<readdy::TimeStep> time;
    std::vector<std::vector<TopologyRecord>> result;

    auto f = h5rd::File::open(filename, h5rd::File::Flag::READ_ONLY);
    if(f->exists(groupName)) {
        auto group = f->getSubgroup(groupName);

        std::size_t nFrames{0};
        // limits
        std::vector<std::size_t> limitsParticles;
        std::vector<std::size_t> limitsEdges;
        {
            if (stride > 1) {
                group.read("limitsParticles", limitsParticles, {stride, 1});
            } else {
                group.read("limitsParticles", limitsParticles);
            }
            if (stride > 1) {
                group.read("limitsEdges", limitsEdges, {stride, 1});
            } else {
                group.read("limitsEdges", limitsEdges);
            }

            if (limitsParticles.size() != limitsEdges.size()) {
                throw std::logic_error(fmt::format("readTopologies: Incompatible limit sizes, "
                                                   "limitsParticles.size={}, limitsEdges.size={}",
                                                   limitsParticles.size(), limitsEdges.size()));
            }

            nFrames = limitsParticles.size() / 2;
            from = std::min(nFrames, from);
            to = std::min(nFrames, to);

            if (from == to) {
                throw std::invalid_argument(fmt::format("readTopologies: not enough frames to cover range ({}, {}]",
                                                        from, to));
            } else {
                limitsParticles = std::vector<std::size_t>(limitsParticles.begin() + 2 * from,
                                                           limitsParticles.begin() + 2 * to);
                limitsEdges = std::vector<std::size_t>(limitsEdges.begin() + 2 * from,
                                                       limitsEdges.begin() + 2 * to);
            }
        }

        from = std::min(nFrames, from);
        to = std::min(nFrames, to);

        // time
        if (stride > 1) {
            group.readSelection("time", time, {from}, {stride}, {to - from});
        } else {
            group.readSelection("time", time, {from}, {1}, {to - from});
        }

        if (limitsParticles.size() % 2 != 0 || limitsEdges.size() % 2 != 0) {
            throw std::logic_error(fmt::format(
                    "limitsParticles size was {} and limitsEdges size was {}, they should be divisible by 2; from={}, to={}",
                    limitsParticles.size(), limitsEdges.size(), from, to)
            );
        }
        // now check that nFrames(particles) == nFrames(edges) == nFrames(time)...
        if (to - from != limitsEdges.size() / 2 || (to - from) != time.size()) {
            throw std::logic_error(fmt::format(
                    "(to-from) should be equal to limitsEdges/2 and equal to the number of time steps in the recording, "
                    "but was: (to-from) = {}, limitsEdges/2 = {}, nTimeSteps={}; from={}, to={}",
                    to - from, limitsEdges.size() / 2, time.size(), from, to
            ));
        }

        std::vector<std::vector<readdy::TopologyTypeId>> types;
        group.readVLENSelection("types", types, {from}, {stride}, {to - from});

        std::size_t ix = 0;
        for (std::size_t frame = from; frame < to; ++frame, ++ix) {
            result.emplace_back();
            // this frame's records
            auto &records = result.back();

            const auto &particlesLimitBegin = limitsParticles.at(2 * ix);
            const auto &particlesLimitEnd = limitsParticles.at(2 * ix + 1);
            // since the edges are flattened, we actually have to multiply this by 2
            const auto &edgesLimitBegin = limitsEdges.at(2 * ix);
            const auto &edgesLimitEnd = limitsEdges.at(2 * ix + 1);

            std::vector<std::size_t> flatParticles;
            group.readSelection("particles", flatParticles, {particlesLimitBegin}, {stride},
                                {particlesLimitEnd - particlesLimitBegin});
            std::vector<std::size_t> flatEdges;
            //readdy::log::critical("edges {} - {}", edgesLimitBegin, edgesLimitEnd);
            group.readSelection("edges", flatEdges, {edgesLimitBegin, 0}, {stride, 1},
                                {edgesLimitEnd - edgesLimitBegin, 2});

            const auto &currentTypes = types.at(ix);
            auto typesIt = currentTypes.begin();
            for (auto particlesIt = flatParticles.begin();
                 particlesIt != flatParticles.end(); ++particlesIt) {

                records.emplace_back();
                auto nParticles = *particlesIt;
                for (std::size_t i = 0; i < nParticles; ++i) {
                    ++particlesIt;
                    records.back().particleIndices.push_back(*particlesIt);
                }
                records.back().type = *typesIt;
                ++typesIt;
            }

            if (currentTypes.size() != records.size()) {
                throw std::logic_error(fmt::format("for frame {} had {} topology types but {} topologies",
                                                   frame, currentTypes.size(), records.size()));
            }

            std::size_t recordIx = 0;
            for (auto edgesIt = flatEdges.begin();
                 edgesIt != flatEdges.end(); ++recordIx) {
                auto &currentRecord = records.at(recordIx);

                auto nEdges = *edgesIt;
                edgesIt += 2;

                for (std::size_t i = 0; i < nEdges; ++i) {
                    currentRecord.edges.emplace_back(*edgesIt, *(edgesIt + 1));
                    edgesIt += 2;
                }

            }

        }

        if (time.size() != result.size()) {
            throw std::logic_error(fmt::format("readTopologies: size mismatch, time size is {}, result size is {}",
                                               time.size(), result.size()));
        }
        return std::make_tuple(time, result);
    }
    return std::make_tuple(time, result);
}

std::vector<std::vector<TrajectoryParticle>> read_trajectory(const std::string &filename, const std::string &name) {
    readdy::io::BloscFilter bloscFilter;
    bloscFilter.registerFilter();

    auto f = h5rd::File::open(filename, h5rd::File::Flag::READ_ONLY);

    auto particleInfoH5Type = readdy::model::ioutils::getParticleTypeInfoType(f->ref());

    // get particle types from config
    std::vector<readdy::model::ioutils::ParticleTypeInfo> types;
    {
        auto config = f->getSubgroup("readdy/config");
        config.read("particle_types", types, &std::get<0>(particleInfoH5Type), &std::get<1>(particleInfoH5Type));
    }
    std::unordered_map<std::size_t, std::string> typeMapping;
    for (const auto &type : types) {
        typeMapping[type.type_id] = std::string(type.name);
    }

    auto traj = f->getSubgroup("readdy/trajectory/" + name);

    // limits
    std::vector<std::size_t> limits;
    traj.read("limits", limits);

    // time
    std::vector<readdy::TimeStep> time;
    traj.read("time", time);

    // records
    std::vector<readdy::model::observables::TrajectoryEntry> entries;
    auto trajectoryEntryTypes = readdy::model::observables::util::getTrajectoryEntryTypes(f->ref());
    traj.read("records", entries, &std::get<0>(trajectoryEntryTypes), &std::get<1>(trajectoryEntryTypes));

    auto n_frames = limits.size() / 2;
    readdy::log::debug("got n frames: {}", n_frames);

    std::vector<std::vector<TrajectoryParticle>> result;
    result.reserve(n_frames);


    auto timeIt = time.begin();
    for (std::size_t frame = 0; frame < limits.size(); frame += 2, ++timeIt) {
        auto begin = limits[frame];
        auto end = limits[frame + 1];
        result.emplace_back();
        auto &currentFrame = result.back();
        currentFrame.reserve(end - begin);

        for (auto it = entries.begin() + begin; it != entries.begin() + end; ++it) {
            currentFrame.emplace_back(typeMapping[it->typeId],
                                      readdy::model::particleflavor::particleFlavorToString(it->flavor),
                                      it->pos.data, it->id, *timeIt);
        }
    }

    return std::move(result);
}

void exportUtils(py::module &m) {
    using namespace pybind11::literals;
    py::class_<TrajectoryParticle>(m, "TrajectoryParticle")
            .def_property_readonly("type", [](const TrajectoryParticle &self) { return self.type; }, R"docs(
                Returns the type of the particle.

                :return: type of the particle
            )docs")
            .def_property_readonly("flavor", [](const TrajectoryParticle &self) { return self.flavor; }, R"docs(
                Returns the flavor of the particle (NORMAL or TOPOLOGY).

                :return: flavor of the particle
            )docs")
            .def_property_readonly("position", [](const TrajectoryParticle &self) { return self.position; }, R"docs(
                Returns the position of the particle as array of length 3.

                :return: position of the particle
            )docs")
            .def_property_readonly("id", [](const TrajectoryParticle &self) { return self.id; }, R"docs(
                Returns the id of the particle.

                :return: id of the particle
            )docs")
            .def_property_readonly("t", [](const TrajectoryParticle &self) { return self.t; }, R"docs(
                Returns the current simulation time.

                :return: the simulation time
            )docs")
            .def("__repr__", [](const TrajectoryParticle &p) {
                return repr(p);
            })
            .def("__str__", [](const TrajectoryParticle &p) {
                return repr(p);
            });
    py::class_<rpy::ReadableReactionRecord>(m, "ReactionRecord")
            .def_property_readonly("type", [](const rpy::ReadableReactionRecord &self) { return self.type; }, R"docs(
                Returns the type of reaction that occurred. One of conversion, fission, fusion, enzymatic, decay.

                :return: the type of reaction
            )docs")
            .def_property_readonly("reaction_label", [](const rpy::ReadableReactionRecord &self) {
                return self.reaction_label;
            }, R"docs(
                Returns the label of the specific reaction as defined in the reaction diffusion system.

                :return: the label
            )docs")
            .def_property_readonly("educts", [](const rpy::ReadableReactionRecord &self) {
                return self.educts;
            }, R"docs(
                Returns the particle IDs of the the educts of this reaction.

                :return: the IDs
            )docs")
            .def_property_readonly("products", [](const rpy::ReadableReactionRecord &self) {
                return self.products;
            }, R"docs(
                Returns the product IDs of the products of this reaction.

                :return: the IDs
            )docs")
            .def_property_readonly("position", [](const rpy::ReadableReactionRecord &self) {
                return self.where;
            }, R"docs(
                Returns the position of the reaction event.

                :return: the position
            )docs")
            .def("__repr__", [](const rpy::ReadableReactionRecord &self) {
                return repr(self);
            })
            .def("__str__", [](const rpy::ReadableReactionRecord &self) {
                return repr(self);
            });
    m.def("convert_xyz", &convert_xyz, "h5_file_name"_a, "traj_data_set_name"_a, "xyz_out_file_name"_a,
          "generate_tcl"_a = true, "tcl_with_grid"_a = false, "radii"_a = radiusmap{},
          "color_ids"_a = std::unordered_map<std::string, unsigned int>{},
          "box_size"_a = std::array<readdy::scalar, 3>{{0.,0.,0.}});
    m.def("convert_readdyviewer", &convert_readdy_viewer, "h5_file_name"_a, "traj_data_set_name"_a,
          "begin"_a = 0, "end"_a = std::numeric_limits<int>::max(), "stride"_a = 1);
    m.def("read_trajectory", &read_trajectory, "filename"_a, "name"_a);
    m.def("read_topologies_observable", &readTopologies, "filename"_a, "groupname"_a,
          "begin"_a = 0, "end"_a = std::numeric_limits<int>::max(), "stride"_a = 1);
    m.def("read_reaction_observable", &read_reactions_obs, "filename"_a, "name"_a);
    py::add_ostream_redirect(m, "ostream_redirect");
}
