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
#pragma once

#include "../DataSet.h"

#include <sstream>
#include <iterator>
#include <atomic>
#include <H5Ppublic.h>
#include <readdy/common/logging.h>
#include <readdy/io/DataSetType.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)
template<typename T, bool VLEN, int compression>
inline void DataSet<T, VLEN, compression>::close() {
    if (memorySpace >= 0 && H5Sclose(memorySpace) < 0) {
        throw std::runtime_error("error on closing memory space!");
    } else {
        memorySpace = -1;
    }
    if (dataSetHandle >= 0 && H5Dclose(dataSetHandle) < 0) {
        throw std::runtime_error("error on closing data set!");
    } else {
        dataSetHandle = -1;
    }
}

template<typename T, bool VLEN, int compression>
inline DataSet<T, VLEN, compression>::~DataSet() {
    close();
}

template<typename T, bool VLEN, int compression>
inline DataSet<T, VLEN, compression>::DataSet(const std::string &name, const Group &group,
                                              const std::vector<h5::dims_t> &chunkSize,
                                              const std::vector<h5::dims_t> &maxDims)
        : DataSet(name, group, chunkSize, maxDims, STDDataSetType<T>(), NativeDataSetType<T>()) {};

template<typename T, bool VLEN, int compression>
inline DataSet<T, VLEN, compression>::DataSet(const std::string &name, const Group &group,
                                              const std::vector<h5::dims_t> &chunkSize,
                                              const std::vector<h5::dims_t> &maxDims, DataSetType memoryType,
                                              DataSetType fileType)
        : maxDims(maxDims), group(group), memoryType(memoryType), fileType(fileType) {
    {
        std::stringstream result;
        std::copy(maxDims.begin(), maxDims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("creating data set (vlen={}) with maxDims={}", VLEN, result.str());
    }
    if (compression == DataSetCompression::blosc) {
        blosc_compression::initialize();
    }
    if (VLEN) {
        log::trace("making the types {} and {} vlen", memoryType.tid, fileType.tid);
        DataSetType vlenMemoryType;
        vlenMemoryType.tid = std::make_shared<DataTypeHandle>(H5Tvlen_create(memoryType.tid->tid));
        this->memoryType = vlenMemoryType;
        DataSetType vlenFileType;
        vlenFileType.tid = std::make_shared<DataTypeHandle>(H5Tvlen_create(fileType.tid->tid));
        this->fileType = vlenFileType;
    }
    // validate and find extension dim
    {
        auto unlimited_it = std::find(maxDims.begin(), maxDims.end(), h5::UNLIMITED_DIMS);
        bool containsUnlimited = unlimited_it != maxDims.end();
        if (!containsUnlimited) {
            throw std::runtime_error("needs to contain unlimited_dims in some dimension to be extensible");
        }
        extensionDim = static_cast<h5::dims_t>(std::distance(maxDims.begin(), unlimited_it));
        log::trace("found extension dim {}", extensionDim);
    }
    {
        // set up empty data set
        std::vector<h5::dims_t> dims(maxDims.begin(), maxDims.end());
        dims[extensionDim] = 0;
        auto fileSpace = H5Screate_simple(static_cast<int>(maxDims.size()), dims.data(), maxDims.data());
        auto plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_chunk(plist, static_cast<int>(chunkSize.size()), chunkSize.data());
        if (compression == DataSetCompression::blosc) {
            blosc_compression::activate(plist);
        }
        dataSetHandle = H5Dcreate(group.handle, name.c_str(), this->fileType.tid->tid, fileSpace, H5P_DEFAULT, plist,
                                  H5P_DEFAULT);
        if (dataSetHandle < 0) {
            log::error("Error on creating data set {}", dataSetHandle);
            H5Eprint(H5Eget_current_stack(), stderr);
            throw std::runtime_error("Error on creating data set " + std::to_string(dataSetHandle));
        }
        memorySpace = -1;
        if (plist >= 0 && H5Pclose(plist) < 0) {
            log::error("Error on closing plist {}", plist);
            H5Eprint(H5Eget_current_stack(), stderr);
            throw std::runtime_error("Error on closing plist " + std::to_string(plist));
        }
        if (fileSpace >= 0 && H5Sclose(fileSpace) < 0) {
            log::error("Error on closing file space {}", fileSpace);
            H5Eprint(H5Eget_current_stack(), stderr);
            throw std::runtime_error("(constructor) Error on closing file space " + std::to_string(fileSpace));
        }
    }
}

template<typename T, bool VLEN, int compression>
template<bool no_vlen>
inline void
DataSet<T, VLEN, compression>::append(std::vector<T> &data, typename std::enable_if<no_vlen, bool>::type *) {
    if (!data.empty()) append({1, data.size()}, data.data());
}

template<typename T, bool VLEN, int compression>
template<bool vlen>
inline void
DataSet<T, VLEN, compression>::append(std::vector<std::vector<T>> &data, typename std::enable_if<vlen, bool>::type *) {
    if (!data.empty()) append({data.size()}, data.data());
}

template<typename T, bool VLEN, int compression>
template<bool enabled>
inline void DataSet<T, VLEN, compression>::append(const std::vector<h5::dims_t> &dims, std::vector<T> *const data,
                                                  typename std::enable_if<enabled, bool>::type *) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("appending to vlen data set with data size = ({})", result.str());
    }
    if (dims.size() != maxDims.size()) {
        log::error("Tried to append data with ndims={} to set with ndims={}", dims.size(), maxDims.size());
        throw std::runtime_error("tried to append data with wrong dimensionality!");
    }
    if (memorySpace == -1) {
        memorySpace = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
    } else {
        H5Sset_extent_simple(memorySpace, static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    std::vector<h5::dims_t> currentExtent;
    {
        currentExtent.resize(dims.size());
        auto fileSpace = H5Dget_space(dataSetHandle);
        H5Sget_simple_extent_dims(fileSpace, currentExtent.data(), nullptr);
        if (fileSpace >= 0 && H5Sclose(fileSpace) < 0) {
            log::error("error on closing file space! {}", fileSpace);
            H5Eprint(H5Eget_current_stack(), stderr);
            throw std::runtime_error("error on closing file space!");
        }
    }
    std::vector<h5::dims_t> offset;
    {
        offset.resize(dims.size());
        offset[extensionDim] = currentExtent[extensionDim];
    }
    std::vector<h5::dims_t> newExtent(currentExtent);
    {
        newExtent[extensionDim] += dims[extensionDim];
        H5Dset_extent(dataSetHandle, newExtent.data());
        {
            log::trace("setting new extent to:");
            std::stringstream result;
            std::copy(newExtent.begin(), newExtent.end(), std::ostream_iterator<int>(result, ", "));
            log::trace("    size = {}", result.str());
        }
    }
    log::trace("selecting hyperslab with:");
    {
        std::stringstream result;
        std::copy(offset.begin(), offset.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("    current extent = {}", result.str());
    }
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("    size = {}", result.str());
    }
    std::vector<hvl_t> traj;
    {
        const auto n = dims[extensionDim];
        log::trace("appending n={} vlen entries", n);
        traj.reserve(n);
        for (auto i = 0; i < n; ++i) {
            auto val = &data[i];
            hvl_t entry;
            entry.len = val->size();
            entry.p = val->data();
            traj.push_back(std::move(entry));
        }
    }
    {
        auto fileSpace = H5Dget_space(dataSetHandle);
        H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
        H5Dwrite(dataSetHandle, memoryType.tid->tid, memorySpace, fileSpace, H5P_DEFAULT, traj.data());
        H5Sclose(fileSpace);
    }
};

template<typename T, bool VLEN, int compression>
template<bool enabled>
inline void DataSet<T, VLEN, compression>::append(const std::vector<h5::dims_t> &dims, const T *const data,
                                                  typename std::enable_if<enabled, bool>::type *) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("appending to regular data set with data size = ({})", result.str());
    }
    if (dims.size() != maxDims.size()) {
        log::error("Tried to append data with ndims={} to set with ndims={}", dims.size(), maxDims.size());
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
            auto fileSpace = H5Dget_space(dataSetHandle);
            H5Sget_simple_extent_dims(fileSpace, currentExtent.data(), nullptr);
            if (fileSpace >= 0 && H5Sclose(fileSpace) < 0) {
                log::error("error on closing file space! {}", fileSpace);
                H5Eprint(H5Eget_current_stack(), stderr);
                throw std::runtime_error("error on closing file space!");
            }
        }
    }
    offset[extensionDim] = currentExtent[extensionDim];
    {
        std::vector<h5::dims_t> newExtent(currentExtent);
        newExtent[extensionDim] += dims[extensionDim];
        H5Dset_extent(dataSetHandle, newExtent.data());
    }
    log::trace("selecting hyperslab with:");
    {
        std::stringstream result;
        std::copy(offset.begin(), offset.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("    current extent = {}", result.str());
    }
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("    size = {}", result.str());
    }
    auto fileSpace = H5Dget_space(dataSetHandle);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
    if (H5Dwrite(dataSetHandle, memoryType.tid->tid, memorySpace, fileSpace, H5P_DEFAULT, data) < 0) {
        log::error("Error with data set {}", dataSetHandle);
        H5Eprint(H5Eget_current_stack(), stderr);
    }
    H5Sclose(fileSpace);
}

template<typename T, bool VLEN, int compression>
inline void DataSet<T, VLEN, compression>::flush() {
    if (dataSetHandle >= 0 && H5Fflush(dataSetHandle, H5F_SCOPE_LOCAL) < 0) {
        throw std::runtime_error("error when flushing HDF5 data set with handle " + std::to_string(dataSetHandle));
    }
}

NAMESPACE_END(io)
NAMESPACE_END(readdy)
