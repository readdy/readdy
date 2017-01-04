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
 * @file DataSet_bits.h
 * @brief << brief description >>
 * @author clonker
 * @date 04/01/2017
 * @copyright GNU Lesser General Public License v3.0
 */
#ifndef READDY_MAIN_DATASET_BITS_H
#define READDY_MAIN_DATASET_BITS_H

#include "../DataSet.h"

#include <H5Ppublic.h>
#include <readdy/common/logging.h>
#include <readdy/io/DataSetType.h>

namespace readdy {
namespace io {
template<typename T, bool VLEN>
inline void DataSet<T, VLEN>::close() {
    if (memorySpace >= 0 && H5Sclose(memorySpace) < 0) {
        throw std::runtime_error("error on closing memory space!");
    } else {
        memorySpace = -1;
    }
    if (handle >= 0 && H5Dclose(handle) < 0) {
        throw std::runtime_error("error on closing data set!");
    } else {
        handle = -1;
    }
}

template<typename T, bool VLEN>
inline DataSet<T, VLEN>::~DataSet() {
    close();
}

template<typename T, bool VLEN>
inline DataSet<T, VLEN>::DataSet(const std::string &name, const Group &group, const std::vector<h5::dims_t> &chunkSize,
                                 const std::vector<h5::dims_t> &maxDims)
        : DataSet(name, group, chunkSize, maxDims, STDDataSetType<T>(), NativeDataSetType<T>()) { };

template<typename T, bool VLEN>
inline DataSet<T, VLEN>::DataSet(const std::string &name, const Group &group, const std::vector<h5::dims_t> &chunkSize,
                                 const std::vector<h5::dims_t> &maxDims, DataSetType memoryType, DataSetType fileType)
        : maxDims(maxDims), group(group), memoryType(memoryType), fileType(fileType) {
    // validate and find extension dim
    {
        auto unlimited_it = std::find(maxDims.begin(), maxDims.end(), h5::UNLIMITED_DIMS);
        bool containsUnlimited = unlimited_it != maxDims.end();
        if (!containsUnlimited) {
            throw std::runtime_error("needs to contain unlimited_dims in some dimension to be extensible");
        }
        extensionDim = static_cast<h5::dims_t>(std::distance(maxDims.begin(), unlimited_it));
        log::console()->trace("found extension dim {}", extensionDim);
    }
    {
        // set up empty data set
        std::vector<h5::dims_t> dims(maxDims.begin(), maxDims.end());
        dims[extensionDim] = 0;
        auto fileSpace = H5Screate_simple(static_cast<int>(maxDims.size()), dims.data(), maxDims.data());
        auto plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_chunk(plist, static_cast<int>(chunkSize.size()), chunkSize.data());
        handle = H5Dcreate(group.handle, name.c_str(), fileType.tid, fileSpace, H5P_DEFAULT, plist, H5P_DEFAULT);
        memorySpace = -1;
        H5Pclose(plist);
        H5Pclose(fileSpace);
    }
}

template<typename T, bool VLEN>
template<bool enabled>
inline void DataSet<T, VLEN>::append(std::vector<T> &data, typename std::enable_if<enabled, bool>::type*) {
    append({1, data.size()}, data.data());
}

template<typename T, bool VLEN>
template<bool enabled>
inline void DataSet<T, VLEN>::append(std::vector<std::vector<T>> &data, typename std::enable_if<enabled, bool>::type*) {
    append({data.size()}, data.data());
}

template<typename T, bool>
struct writer { };


template<typename T>
struct writer<T, true> {
    static void
    write(h5::handle_t handle, h5::data_set_type_t memType, h5::handle_t memSpace, hid_t fSpace, const std::vector<T> *const data,
          h5::dims_t extensionDim, const std::vector<h5::dims_t> &dims) {
        std::vector<hvl_t> traj;
        {
            const auto n = dims[extensionDim];
            traj.reserve(n);
            for (auto i = 0; i < n; ++i) {
                auto val = data[i];
                hvl_t entry;
                entry.len = val.size();
                entry.p = val.data();
                traj.push_back(std::move(entry));
            }
        }
        H5Dwrite(handle, memType, memSpace, fSpace, H5P_DEFAULT, traj.data());
    };
};

template<typename T>
struct writer<T, false> {
    static void
    write(h5::handle_t handle, h5::data_set_type_t memType, h5::handle_t memSpace, hid_t fSpace, const T *const data,
          h5::dims_t extensionDim, const std::vector<h5::dims_t> &dims) {
        H5Dwrite(handle, memType, memSpace, fSpace, H5P_DEFAULT, data);
    };
};

template<typename T, bool VLEN>
template<bool enabled>
inline void DataSet<T, VLEN>::append(const std::vector<h5::dims_t> &dims, const std::vector<T> *const data, typename std::enable_if<enabled, bool>::type*) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::console()->trace("appending to vlen data set with data size = ({})", result.str());
    }
    if (dims.size() != maxDims.size()) {
        log::console()->error("Tried to append data with ndims={} to set with ndims={}", dims.size(), maxDims.size());
        throw std::runtime_error("tried to append data with wrong dimensionality!");
    }
    if (memorySpace == -1) {
        memorySpace = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
    } else {
        H5Sset_extent_simple(memorySpace, static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    std::vector<h5::dims_t> currentExtent;
    std::vector<h5::dims_t> offset;
    offset.resize(dims.size());
    {
        currentExtent.resize(dims.size());
        {
            auto fileSpace = H5Dget_space(handle);
            H5Sget_simple_extent_dims(fileSpace, currentExtent.data(), nullptr);
            H5Dclose(fileSpace);
        }
    }
    offset[extensionDim] = currentExtent[extensionDim];
    {
        std::vector<h5::dims_t> newExtent(currentExtent);
        newExtent[extensionDim] += dims[extensionDim];
        H5Dset_extent(handle, newExtent.data());
    }
    auto fileSpace = H5Dget_space(handle);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
    writer<T, VLEN>::write(handle, memoryType.tid, memorySpace, fileSpace, data, extensionDim, dims);
    H5Sclose(fileSpace);
};

template<typename T, bool VLEN>
template<bool enabled>
inline void DataSet<T, VLEN>::append(const std::vector<h5::dims_t> &dims, const T *const data, typename std::enable_if<enabled, bool>::type*) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::console()->trace("appending to regular data set with data size = ({})", result.str());
    }
    if (dims.size() != maxDims.size()) {
        log::console()->error("Tried to append data with ndims={} to set with ndims={}", dims.size(), maxDims.size());
        throw std::runtime_error("tried to append data with wrong dimensionality!");
    }
    if (memorySpace == -1) {
        memorySpace = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
    } else {
        H5Sset_extent_simple(memorySpace, static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    std::vector<h5::dims_t> currentExtent;
    std::vector<h5::dims_t> offset;
    offset.resize(dims.size());
    {
        currentExtent.resize(dims.size());
        {
            auto fileSpace = H5Dget_space(handle);
            H5Sget_simple_extent_dims(fileSpace, currentExtent.data(), nullptr);
            H5Dclose(fileSpace);
        }
    }
    offset[extensionDim] = currentExtent[extensionDim];
    {
        std::vector<h5::dims_t> newExtent(currentExtent);
        newExtent[extensionDim] += dims[extensionDim];
        H5Dset_extent(handle, newExtent.data());
    }
    auto fileSpace = H5Dget_space(handle);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
    //H5Dwrite(handle, memoryType.tid, memorySpace, fileSpace, H5P_DEFAULT, data);
    writer<T, VLEN>::write(handle, memoryType.tid, memorySpace, fileSpace, data, extensionDim, dims);
    H5Sclose(fileSpace);
}

}
}
#endif //READDY_MAIN_DATASET_BITS_H
