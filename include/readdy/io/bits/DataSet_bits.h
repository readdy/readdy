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

#include <sstream>
#include <iterator>
#include <atomic>

#include <H5Ppublic.h>
#include <readdy/common/logging.h>

#include "../DataSet.h"
#include "../DataSetType.h"
#include "../PropertyList.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

inline DataSpace DataSet::getFileSpace() const {
    auto _hid = H5Dget_space(hid());
    if (_hid < 0) {
        log::error("Failed to get data set space!");
        H5Eprint(H5Eget_current_stack(), stderr);
    }
    return DataSpace(_hid);
}

inline DataSet::~DataSet() = default;

inline DataSet::DataSet(h5::h5_handle handle, const DataSetType &memoryType, const DataSetType &fileType)
        : Object(std::make_shared<DataSetHandle>(handle)), memoryType(memoryType), fileType(fileType),
          memorySpace(-1) {}

template<typename T>
inline void DataSet::append(std::vector<T> &data) {
    if (!data.empty()) append({1, data.size()}, data.data());
}

template<typename T>
inline void VLENDataSet::append(std::vector<std::vector<T>> &data) {
    if (!data.empty()) append({data.size()}, data.data());
}

template<typename T>
inline void VLENDataSet::append(const std::vector<h5::h5_dims> &dims, std::vector<T> *const data) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("appending to vlen data set with data size = ({})", result.str());
    }
    {
        auto fs = getFileSpace();
        if (dims.size() != fs.ndim()) {
            log::error("Tried to append data with ndims={} to set with ndims={}", dims.size(), fs.ndim());
            throw std::runtime_error("tried to append data with wrong dimensionality!");
        }
    }
    if (memorySpace.hid() < 0) {
        memorySpace = DataSpace(dims);
    } else {
        H5Sset_extent_simple(memorySpace.hid(), static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    std::vector<h5::h5_dims> currentExtent;
    {
        auto fs = getFileSpace();
        currentExtent = fs.dims();
    }
    std::vector<h5::h5_dims> offset;
    {
        offset.resize(dims.size());
        offset[_extensionDim] = currentExtent[_extensionDim];
    }
    std::vector<h5::h5_dims> newExtent(currentExtent.begin(), currentExtent.end());
    {
        newExtent[_extensionDim] += dims[_extensionDim];
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
        const auto n = dims[_extensionDim];
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
        H5Dwrite(hid(), memoryType.hid(), memorySpace.hid(), fileSpace.hid(), H5P_DEFAULT, traj.data());
    }
}

template<typename T>
inline void DataSet::append(const std::vector<h5::h5_dims> &dims, const T *const data) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("appending to regular data set with data size = ({})", result.str());
    }
    if (dims.size() != getFileSpace().ndim()) {
        log::error("Tried to append data with ndims={} to set with ndims={}", dims.size(), getFileSpace().ndim());
        throw std::runtime_error("tried to append data with wrong dimensionality!");
    }
    if (memorySpace.hid() < 0) {
        memorySpace = DataSpace(dims);
    } else {
        H5Sset_extent_simple(memorySpace.hid(), static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    std::vector<h5::h5_dims> currentExtent;
    std::vector<h5::h5_dims> offset;
    offset.resize(dims.size());
    {
        currentExtent = getFileSpace().dims();
    }
    offset[_extensionDim] = currentExtent[_extensionDim];
    {
        std::vector<h5::h5_dims> newExtent(currentExtent);
        newExtent[_extensionDim] += dims[_extensionDim];
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
    if (H5Dwrite(hid(), memoryType.hid(), memorySpace.hid(), fileSpace.hid(), H5P_DEFAULT, data) < 0) {
        log::error("Error with data set {}", hid());
        H5Eprint(H5Eget_current_stack(), stderr);
    }
}

inline void DataSet::flush() {
    if (hid() >= 0 && H5Fflush(hid(), H5F_SCOPE_LOCAL) < 0) {
        throw std::runtime_error("error when flushing HDF5 data set with handle " + std::to_string(hid()));
    }
}

inline VLENDataSet::VLENDataSet(h5::h5_handle handle, const DataSetType &memoryType, const DataSetType &fileType)
        : Object(std::make_shared<DataSetHandle>(handle)), memoryType(memoryType), fileType(fileType),
          memorySpace(-1) {
    if(handle < 0) {
        log::critical("tried to create a data set with negative handle!");
    }
}

NAMESPACE_END(io)
NAMESPACE_END(readdy)
