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
 * @file TestIO.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 25.05.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/io/TrajectoryEntry.h>
#include <readdy/model/observables/io/Types.h>
#include <fstream>
#include "gtest/gtest.h"
#include "readdy/readdy.h"

namespace {
struct RAIIFS {
    explicit RAIIFS(const std::string &fname) : fs() {
        fs.open(fname, std::fstream::out | std::fstream::trunc);
    }

    RAIIFS(const RAIIFS &) = delete;

    RAIIFS &operator=(const RAIIFS &) = delete;

    RAIIFS(RAIIFS &&) = delete;

    RAIIFS &operator=(RAIIFS &&) = delete;

    ~RAIIFS() {
        fs.close();
    }

    std::fstream fs;
};

class XYZConverter {
public:

    XYZConverter(std::string h5name, std::string trajName, std::string out)
            : _h5name(std::move(h5name)), _out(std::move(out)), _trajName(std::move(trajName)) {
        readdy::io::blosc_compression::initialize();
    };

    void convert() {
        readdy::log::debug(R"(converting "{}" to "{}")", _h5name, _out);

        readdy::io::File f(_h5name, readdy::io::File::Action::OPEN, readdy::io::File::Flag::READ_ONLY);
        auto &rootGroup = f.getRootGroup();

        // get particle types from config
        std::vector<readdy::model::ioutils::ParticleTypeInfo> types;
        {
            auto config = rootGroup.subgroup("readdy/config");
            config.read("particle_types", types, readdy::model::ioutils::ParticleTypeInfoMemoryType(),
                        readdy::model::ioutils::ParticleTypeInfoFileType());
        }

        // map from type name to max number of particles in traj
        std::unordered_map<readdy::particle_type_type, std::size_t> maxCounts;
        for (const auto &type : types) {
            maxCounts[type.type_id] = 0;
            readdy::log::debug("got type {} with id {} and D {}", type.name, type.type_id, type.diffusion_constant);
        }

        auto traj = rootGroup.subgroup("readdy/trajectory/" + _trajName);

        // limits
        std::vector<std::size_t> limits;
        traj.read("limits", limits, readdy::io::STDDataSetType<std::size_t>(),
                  readdy::io::NativeDataSetType<std::size_t>());
        // records
        std::vector<readdy::model::observables::TrajectoryEntry> entries;
        readdy::model::observables::util::TrajectoryEntryMemoryType memoryType;
        readdy::model::observables::util::TrajectoryEntryFileType fileType;
        traj.read("records", entries, memoryType, fileType);

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

        //std::vector<readdy::time_step_type> timesteps;
        //traj.read("time", timesteps);

        readdy::log::debug("writing to xyz (n timesteps {})", limits.size() / 2);

        {
            RAIIFS fshandle{_out};
            auto &fs = fshandle.fs;

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

                for(const auto &mappingEntry : typeMapping) {
                    auto nGhosts = maxCounts.at(mappingEntry.first) - currentCounts.at(mappingEntry.second);

                    fs << xyzPerType.at(mappingEntry.second);
                    for(int x = 0; x < nGhosts; ++x) {
                        fs << "type_" + std::to_string(mappingEntry.first) + "\t0\t0\t0\n";
                    }
                }
            }

        }

        readdy::log::debug("converting finished");
    }

private:
    std::string _h5name, _out, _trajName;
};

TEST(TestIO, ReadTrajectory) {

    /*using namespace readdy;
    readdy::io::blosc_compression::initialize();
    io::File f("", io::File::Action::OPEN, io::File::Flag::READ_ONLY);
    auto& rootGroup = f.getRootGroup();
    auto config = rootGroup.subgroup("readdy/config");
    std::vector<readdy::model::ioutils::ParticleTypeInfo> types;
    config.read("particle_types", types, readdy::model::ioutils::ParticleTypeInfoMemoryType(), readdy::model::ioutils::ParticleTypeInfoFileType());

    for(const auto& info : types) {
        readdy::log::warn("foo: {}", info.name);
        readdy::log::warn("foo: {}", info.diffusion_constant);
        readdy::log::warn("foo: {}", info.type_id);
    }*/
    readdy::log::console()->set_level(spdlog::level::debug);


}

}