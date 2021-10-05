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
 * @file DataSet_detail.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include <iostream>
#include <iterator>
#include <utility>

#include "../DataSet.h"
#include "../DataSpace.h"
#include "../DataSetType.h"

inline h5rd::DataSet::~DataSet() {
    try {
        close();
    } catch (const Exception &e) {
        std::cerr << "Unable to close hdf5 data set: " << e.what() << std::endl;
    }
}

inline void h5rd::DataSet::close() {
    auto pf = _parentFile.lock();
    if (pf) {
        if (!pf->closed() && valid() && H5Dclose(id()) < 0) {
            throw Exception("Error on closing HDF5 data set");
        }
    }
}

inline h5rd::DataSet::DataSet(ParentFileRef parentFile, const DataSetType &memoryType, const DataSetType &fileType)
        : SubObject(std::move(parentFile)), _memoryType(memoryType), _fileType(fileType) {}

inline h5rd::dimension &h5rd::DataSet::extensionDim() {
    return _extensionDim;
}

inline const h5rd::dimension &h5rd::DataSet::extensionDim() const {
    return _extensionDim;
}

inline std::shared_ptr<h5rd::DataSpace> h5rd::DataSet::getFileSpace() const {
    auto _hid = H5Dget_space(id());
    if (_hid < 0) {
        throw Exception("Failed to get file space for data set!");
    }
    return std::make_shared<h5rd::DataSpace>(_parentFile, _hid);
}

inline void h5rd::DataSet::flush() {
    if (valid() && H5Fflush(id(), H5F_SCOPE_LOCAL) < 0) {
        throw Exception("error when flushing HDF5 data set with handle " + std::to_string(id()));
    }
}

template<typename T>
inline void h5rd::DataSet::append(std::vector<T> &data) {
    if (!data.empty()) append({1, data.size()}, data.data());
}

template<typename T>
inline void h5rd::DataSet::append(const h5rd::dimensions &dims, const T *data) {
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        // todo log::trace("appending to regular data set with data size = ({})", result.str());
    }
    if (dims.size() != getFileSpace()->ndim()) {
        // todo log::error("Tried to append data with ndims={} to set with ndims={}", dims.size(), getFileSpace().ndim());
        throw std::invalid_argument("tried to append data with wrong dimensionality!");
    }
    if (!_memorySpace) {
        _memorySpace = std::make_unique<DataSpace>(_parentFile, dims);
    } else {
        H5Sset_extent_simple(_memorySpace->id(), static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    dimensions currentExtent;
    dimensions offset;
    offset.resize(dims.size());
    {
        currentExtent = getFileSpace()->dims();
    }
    offset[_extensionDim] = currentExtent[_extensionDim];
    {
        dimensions newExtent(currentExtent);
        newExtent[_extensionDim] += dims[_extensionDim];
        H5Dset_extent(id(), newExtent.data());
    }
    // todo log::trace("selecting hyperslab with:");
    {
        std::stringstream result;
        std::copy(offset.begin(), offset.end(), std::ostream_iterator<int>(result, ", "));
        // todo log::trace("    current extent = {}", result.str());
    }
    {
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        // todo log::trace("    size = {}", result.str());
    }
    auto fileSpace = getFileSpace();
    H5Sselect_hyperslab(fileSpace->id(), H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
    if (H5Dwrite(id(), _memoryType.id(), _memorySpace->id(), fileSpace->id(), H5P_DEFAULT, data) < 0) {
        //log::error("Error with data set {}", hid());
        //H5Eprint(H5Eget_current_stack(), stderr);
        throw Exception("Error on writing data set " + std::to_string(id()));
    }
}
