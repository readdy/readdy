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
template<typename T, bool VLEN>
inline void DataSet<T, VLEN>::close() {
    memorySpaceRef.reset();
}

template<typename T, bool vlen>
inline DataSpace DataSet<T, vlen>::getFileSpace() const {
    auto _hid = H5Dget_space(hid());
    if (_hid < 0) {
        log::error("Failed to get data set space!");
        H5Eprint(H5Eget_current_stack(), stderr);
    }
    return DataSpace(_hid);
}

template<typename T, bool VLEN>
inline DataSet<T, VLEN>::~DataSet() = default;

template<typename T, bool VLEN>
inline DataSet<T, VLEN>::DataSet(h5::handle_t handle, DataSetType memoryType, DataSetType fileType)
        : Object(std::make_shared<DataSetHandle>(handle)), memoryType(memoryType), fileType(fileType) {}

template<typename T, bool VLEN>
inline DataSet<T, VLEN>::DataSet(const std::string &name, const Group &group,
                                              const std::vector<h5::dims_t> &chunkSize,
                                              const std::vector<h5::dims_t> &maxDims, DataSetCompression compression)
        : DataSet(name, group, chunkSize, maxDims, STDDataSetType<T>(), NativeDataSetType<T>(), compression) {};

template<typename T, bool VLEN>
inline DataSet<T, VLEN>::DataSet(const std::string &name, const Group &group,
                                              const std::vector<h5::dims_t> &chunkSize,
                                              const std::vector<h5::dims_t> &maxDims, DataSetType memoryType,
                                              DataSetType fileType, DataSetCompression compression)
        : Object(std::make_shared<DataSetHandle>(-1)), memoryType(memoryType), fileType(fileType) {
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
        DataSpace fileSpace(dims, maxDims);
        auto plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_chunk(plist, static_cast<int>(chunkSize.size()), chunkSize.data());
        unsigned int cd_values[7];
        if (compression == DataSetCompression::blosc) {
            blosc_compression::activate(plist, cd_values);
        }
        auto hid = H5Dcreate(group.hid(), name.c_str(), this->fileType.tid->tid,
                             fileSpace.hid(), H5P_DEFAULT, plist, H5P_DEFAULT);
        if (hid < 0) {
            log::error("Error on creating data set {}", hid);
            H5Eprint(H5Eget_current_stack(), stderr);
            throw std::runtime_error("Error on creating data set " + std::to_string(hid));
        } else {
            handle->set(hid);
        }
        if (plist >= 0 && H5Pclose(plist) < 0) {
            log::error("Error on closing plist {}", plist);
            H5Eprint(H5Eget_current_stack(), stderr);
            throw std::runtime_error("Error on closing plist " + std::to_string(plist));
        }
    }
}

template<typename T, bool VLEN>
template<bool no_vlen>
inline void
DataSet<T, VLEN>::append(std::vector<T> &data, typename std::enable_if<no_vlen, bool>::type *) {
    if (!data.empty()) append({1, data.size()}, data.data());
}

template<typename T, bool VLEN>
template<bool vlen>
inline void
DataSet<T, VLEN>::append(std::vector<std::vector<T>> &data, typename std::enable_if<vlen, bool>::type *) {
    if (!data.empty()) append({data.size()}, data.data());
}

template<typename T, bool VLEN>
template<bool enabled>
inline void DataSet<T, VLEN>::append(const std::vector<h5::dims_t> &dims, std::vector<T> *const data,
                                                  typename std::enable_if<enabled, bool>::type *) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("appending to vlen data set with data size = ({})", result.str());
    }
    if (dims.size() != getFileSpace().ndim()) {
        log::error("Tried to append data with ndims={} to set with ndims={}", dims.size(), getFileSpace().ndim());
        throw std::runtime_error("tried to append data with wrong dimensionality!");
    }
    if (!memorySpaceRef) {
        memorySpaceRef = std::make_unique<DataSpace>(dims);
    } else {
        H5Sset_extent_simple(memorySpaceRef->hid(), static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    std::vector<h5::dims_t> currentExtent;
    {
        currentExtent = getFileSpace().dims();
    }
    std::vector<h5::dims_t> offset;
    {
        offset.resize(dims.size());
        offset[extensionDim] = currentExtent[extensionDim];
    }
    std::vector<h5::dims_t> newExtent(currentExtent.begin(), currentExtent.end());
    {
        newExtent[extensionDim] += dims[extensionDim];
        H5Dset_extent(hid(), newExtent.data());
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
        auto fileSpace = getFileSpace();
        H5Sselect_hyperslab(fileSpace.hid(), H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
        H5Dwrite(hid(), memoryType.tid->tid, memorySpaceRef->hid(), fileSpace.hid(), H5P_DEFAULT,
                 traj.data());
    }
};

template<typename T, bool VLEN>
template<bool enabled>
inline void DataSet<T, VLEN>::append(const std::vector<h5::dims_t> &dims, const T *const data,
                                                  typename std::enable_if<enabled, bool>::type *) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("appending to regular data set with data size = ({})", result.str());
    }
    if (dims.size() != getFileSpace().ndim()) {
        log::error("Tried to append data with ndims={} to set with ndims={}", dims.size(), getFileSpace().ndim());
        throw std::runtime_error("tried to append data with wrong dimensionality!");
    }
    if (!memorySpaceRef) {
        memorySpaceRef = std::make_unique<DataSpace>(dims);
    } else {
        H5Sset_extent_simple(memorySpaceRef->hid(), static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    std::vector<h5::dims_t> currentExtent;
    std::vector<h5::dims_t> offset;
    offset.resize(dims.size());
    {
        currentExtent = getFileSpace().dims();
    }
    offset[extensionDim] = currentExtent[extensionDim];
    {
        std::vector<h5::dims_t> newExtent(currentExtent);
        newExtent[extensionDim] += dims[extensionDim];
        H5Dset_extent(hid(), newExtent.data());
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
    auto fileSpace = getFileSpace();
    H5Sselect_hyperslab(fileSpace.hid(), H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
    if (H5Dwrite(hid(), memoryType.tid->tid, memorySpaceRef->hid(), fileSpace.hid(), H5P_DEFAULT, data) <
        0) {
        log::error("Error with data set {}", hid());
        H5Eprint(H5Eget_current_stack(), stderr);
    }
}

template<typename T, bool VLEN>
inline void DataSet<T, VLEN>::flush() {
    if (hid() >= 0 && H5Fflush(hid(), H5F_SCOPE_LOCAL) < 0) {
        throw std::runtime_error("error when flushing HDF5 data set with handle " + std::to_string(hid()));
    }
}

NAMESPACE_END(io)
NAMESPACE_END(readdy)
