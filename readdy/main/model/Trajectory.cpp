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
#include <hdf5.h>

namespace io = readdy::io;

using particle_t = readdy::model::Particle;
using vec_t = readdy::model::Vec3;

namespace readdy {
namespace model {
namespace observables {

const std::string Trajectory::TRAJECTORY_GROUP_PATH = "/readdy";
const std::string Trajectory::TRAJECTORY_DATA_SET_NAME = "/readdy/trajectory";

Trajectory::Trajectory(readdy::model::Kernel *const kernel, unsigned int stride, unsigned int fs, io::File &file)
        : base_t(kernel, stride), flushStride(fs), file(file) {
    auto trajectoryGroup = file.createGroup(TRAJECTORY_GROUP_PATH);
    entriesTypeMemory = H5Tvlen_create(getEntryTypeMemory());
    entriesTypeFile = H5Tvlen_create(getEntryTypeFile());
    std::vector<unsigned long long> vec{0};
    auto fileSpace = H5Screate_simple(static_cast<int>(1), vec.data(), &io::Object::UNLIMITED_DIMS);
    auto plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    {
        io::Object::dims_t _fs = flushStride;
        H5Pset_chunk(plist, 1, &_fs);
    }
    dataSetHandle = H5Dcreate(trajectoryGroup.getHandle(), TRAJECTORY_DATA_SET_NAME.c_str(),
                              entriesTypeFile, fileSpace, H5P_DEFAULT, plist, H5P_DEFAULT);
    memorySpace = -1;
    H5Pclose(plist);
    H5Sclose(fileSpace);
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
    file.flush();
    result.clear();
    count = 0;
}

void Trajectory::append(observable_entry_t &result) {
    const hsize_t size = result.size();
    if (size > 0) {
        // resize / create memory space
        {
            if (memorySpace == -1) {
                memorySpace = H5Screate_simple(1, &size, nullptr);
            } else {
                H5Sset_extent_simple(memorySpace, 1, &size, nullptr);
            }
        }
        std::vector<hvl_t> traj;
        {
            traj.reserve(result.size());
            for (auto &val : result) {
                hvl_t entry;
                entry.len = val.size();
                entry.p = val.data();
                traj.push_back(std::move(entry));
            }
        }
        io::Object::dims_t currentExtent;
        {
            auto fileSpace = H5Dget_space(dataSetHandle);
            H5Sget_simple_extent_dims(fileSpace, &currentExtent, nullptr);
            H5Sclose(fileSpace);
        }
        {
            io::Object::dims_t newExtent = currentExtent + result.size();
            H5Dset_extent(dataSetHandle, &newExtent);
        }
        auto fileSpace = H5Dget_space(dataSetHandle);
        H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, &currentExtent, nullptr, &size, nullptr);
        H5Dwrite(dataSetHandle, entriesTypeMemory, memorySpace, fileSpace, H5P_DEFAULT, traj.data());
        H5Sclose(fileSpace);
    }
}

io::Object::data_set_type_t Trajectory::getEntryTypeMemory() {
    static hid_t entryTypeMemory = H5Tcreate(H5T_COMPOUND, sizeof(TrajectoryEntry));
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        // init vec pod
        // init trajectory entry pod
        H5Tinsert(entryTypeMemory, "typeId", HOFFSET(TrajectoryEntry, typeId),
                  io::NativeDataSetType<particle_t::type_type>::tid);
        H5Tinsert(entryTypeMemory, "t", HOFFSET(TrajectoryEntry, t),
                  io::NativeDataSetType<readdy::model::observables::time_step_type>::tid);
        H5Tinsert(entryTypeMemory, "id", HOFFSET(TrajectoryEntry, id), io::NativeDataSetType<particle_t::id_type>::tid);
        H5Tinsert(entryTypeMemory, "px", HOFFSET(TrajectoryEntry, px), io::NativeDataSetType<TrajectoryEntry::pos_t>::tid);
        H5Tinsert(entryTypeMemory, "py", HOFFSET(TrajectoryEntry, py), io::NativeDataSetType<TrajectoryEntry::pos_t>::tid);
        H5Tinsert(entryTypeMemory, "pz", HOFFSET(TrajectoryEntry, pz), io::NativeDataSetType<TrajectoryEntry::pos_t>::tid);
    }
    return entryTypeMemory;
}

io::Object::data_set_type_t Trajectory::getEntryTypeFile() {
    static bool initialized = false;
    std::size_t size = sizeof(particle_t::type_type) + sizeof(readdy::model::observables::time_step_type) +
                       sizeof(particle_t::id_type) + 3 * sizeof(TrajectoryEntry::pos_t);
    static hid_t entryTypeFile = H5Tcreate(H5T_COMPOUND, size);

    if (!initialized) {
        initialized = true;
        // init trajectory pod
        H5Tinsert(entryTypeFile, "typeId", 0, io::STDDataSetType<particle_t::type_type>::tid);
        auto cumsize = sizeof(particle_t::type_type);
        H5Tinsert(entryTypeFile, "t", cumsize, io::STDDataSetType<readdy::model::observables::time_step_type>::tid);
        cumsize += sizeof(readdy::model::observables::time_step_type);
        H5Tinsert(entryTypeFile, "id", cumsize, io::STDDataSetType<particle_t::id_type>::tid);
        cumsize += sizeof(particle_t::id_type);
        H5Tinsert(entryTypeFile, "px", cumsize, io::STDDataSetType<TrajectoryEntry::pos_t>::tid);
        cumsize += sizeof(TrajectoryEntry::pos_t);
        H5Tinsert(entryTypeFile, "py", cumsize, io::STDDataSetType<TrajectoryEntry::pos_t>::tid);
        cumsize += sizeof(TrajectoryEntry::pos_t);
        H5Tinsert(entryTypeFile, "pz", cumsize, io::STDDataSetType<TrajectoryEntry::pos_t>::tid);
    }
    return entryTypeFile;
}

Trajectory::~Trajectory() {
    flush();
    if (memorySpace >= 0 && H5Sclose(memorySpace) < 0) {
        log::console()->error("error on closing memory space!");
    } else {
        memorySpace = -1;
    }
    if (dataSetHandle >= 0 && H5Dclose(dataSetHandle) < 0) {
        log::console()->error("error on closing data set!");
    } else {
        dataSetHandle = -1;
    }
    if (entriesTypeMemory >= 0 && H5Tclose(entriesTypeMemory) < 0) {
        log::console()->error("error on closing entries type memory");
    }
    if (entriesTypeFile >= 0 && H5Tclose(entriesTypeFile) < 0) {
        log::console()->error("error on closing entries type file");
    }

}

}
}
}