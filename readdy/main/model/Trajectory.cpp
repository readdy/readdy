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
 * @file Trajectory.cpp
 * @brief Core library impl for Trajectory.h
 * @author chrisfroe
 * @date 12.12.16
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/io/Trajectory.h>
#include <readdy/io/File.h>
#include <readdy/io/DataSet.h>

namespace io = readdy::io;

using particle_t = readdy::model::Particle;
using vec_t = readdy::model::Vec3;
namespace readdy {
namespace model {
namespace observables {

std::ostream &operator<<(std::ostream &os, const TrajectoryEntry &p) {
    os << "TrajectoryEntry(id=" << p.id << ")";
    return os;
}

class TrajectoryEntryMemoryType : public readdy::io::DataSetType {
public:
    TrajectoryEntryMemoryType() {
        using entry_t = readdy::model::observables::TrajectoryEntry;
        tid = []() -> hid_t {
            hid_t entryTypeMemory = H5Tcreate(H5T_COMPOUND, sizeof(entry_t));
            // data types
            io::NativeDataSetType<particle_t::type_type> particleTypeType{};
            io::NativeDataSetType<readdy::model::observables::time_step_type> timeStepType{};
            io::NativeDataSetType<particle_t::id_type> particleIdType{};
            io::NativeDataSetType<entry_t::pos_t> posType{};
            // init vec pod
            // init trajectory entry pod
            H5Tinsert(entryTypeMemory, "typeId", HOFFSET(entry_t, typeId), particleTypeType.tid);
            H5Tinsert(entryTypeMemory, "t", HOFFSET(entry_t, t), timeStepType.tid);
            H5Tinsert(entryTypeMemory, "id", HOFFSET(entry_t, id), particleIdType.tid);
            H5Tinsert(entryTypeMemory, "px", HOFFSET(entry_t, px), posType.tid);
            H5Tinsert(entryTypeMemory, "py", HOFFSET(entry_t, py), posType.tid);
            H5Tinsert(entryTypeMemory, "pz", HOFFSET(entry_t, pz), posType.tid);
            return entryTypeMemory;
        }();
        log::console()->debug("created trajectory entry memory type {}", tid);
    }
};

class TrajectoryEntryFileType : public readdy::io::DataSetType {
public:
    TrajectoryEntryFileType() {
        using entry_t = readdy::model::observables::TrajectoryEntry;
        tid = []() -> hid_t {
            std::size_t size = sizeof(particle_t::type_type) + sizeof(readdy::model::observables::time_step_type) +
                               sizeof(particle_t::id_type) + 3 * sizeof(entry_t::pos_t);
            hid_t entryTypeFile = H5Tcreate(H5T_COMPOUND, size);
            // data types
            io::STDDataSetType<particle_t::type_type> particleTypeType{};
            io::STDDataSetType<readdy::model::observables::time_step_type> timeStepType{};
            io::STDDataSetType<particle_t::id_type> particleIdType{};
            io::STDDataSetType<entry_t::pos_t> posType{};
            // init trajectory pod
            H5Tinsert(entryTypeFile, "typeId", 0, particleTypeType.tid);
            auto cumsize = sizeof(particle_t::type_type);
            H5Tinsert(entryTypeFile, "t", cumsize, timeStepType.tid);
            cumsize += sizeof(readdy::model::observables::time_step_type);
            H5Tinsert(entryTypeFile, "id", cumsize, particleIdType.tid);
            cumsize += sizeof(particle_t::id_type);
            H5Tinsert(entryTypeFile, "px", cumsize, posType.tid);
            cumsize += sizeof(entry_t::pos_t);
            H5Tinsert(entryTypeFile, "py", cumsize, posType.tid);
            cumsize += sizeof(entry_t::pos_t);
            H5Tinsert(entryTypeFile, "pz", cumsize, posType.tid);
            return entryTypeFile;
        }();
        log::console()->debug("created trajectory entry file type {}", tid);
    }
};


const std::string Trajectory::TRAJECTORY_GROUP_PATH = "/readdy";
const std::string Trajectory::TRAJECTORY_DATA_SET_NAME = "/readdy/trajectory";

struct Trajectory::Impl {
    readdy::io::File &file;
    readdy::io::DataSet<TrajectoryEntry, true> dataSet;

    Impl(io::File &file, unsigned int flushStride) : file(file),
                                                     dataSet{
                                                             TRAJECTORY_DATA_SET_NAME,
                                                             file.createGroup(TRAJECTORY_GROUP_PATH),
                                                             {flushStride}, {readdy::io::h5::UNLIMITED_DIMS},
                                                             TrajectoryEntryMemoryType(),
                                                             TrajectoryEntryFileType()
                                                     } {}
};


Trajectory::Trajectory(readdy::model::Kernel *const kernel, unsigned int stride, unsigned int fs, io::File &file)
        : base_t(kernel, stride), pimpl(std::make_unique<Impl>(file, fs)), flushStride(fs) {
};

void Trajectory::evaluate() {
    {
        std::vector<TrajectoryEntry> entries;
        {
            const auto &currentInput = kernel->getKernelStateModel().getParticles();
            entries.reserve(currentInput.size());
            const auto t = getCurrentTimeStep();
            std::transform(currentInput.begin(), currentInput.end(), std::back_inserter(entries),
                           [t](const particle_t &p) { return TrajectoryEntry(p, t); });
        }
        result.reserve(flushStride > 0 ? flushStride : 1);
        result.push_back(std::move(entries));
    }
    if (flushStride == 0 || count % flushStride == 0) {
        flush();
    }
    ++count;
}

void Trajectory::flush() {
    append(result);
    pimpl->file.flush();
    result.clear();
    count = 0;
}

void Trajectory::append(observable_entry_t &result) {
    pimpl->dataSet.append(result);
}

Trajectory::~Trajectory() {
    flush();
}

void Trajectory::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    throw std::runtime_error("not supported for trajectory");
}

void Trajectory::append() {
    throw std::runtime_error("not supported for trajectory");
}

}
}
}