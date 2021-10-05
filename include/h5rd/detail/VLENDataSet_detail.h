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
 * @file VLENDataSet.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright BSD-3
 */

#pragma once

#include "../VLENDataSet.h"
#include "../DataSpace.h"

namespace h5rd {

inline VLENDataSet::VLENDataSet(ParentFileRef parentFile, const DataSetType &memoryType, const DataSetType &fileType)
        : SubObject(std::move(parentFile)), _memoryType(memoryType), _fileType(fileType) {}

inline std::shared_ptr<DataSpace> VLENDataSet::getFileSpace() const {
    auto _hid = H5Dget_space(id());
    if (_hid < 0) {
        throw Exception("Failed to get file space for vlen data set!");
    }
    return std::make_shared<h5rd::DataSpace>(_parentFile, _hid);
}

inline dimension &VLENDataSet::extensionDim() {
    return _extensionDim;
}

inline const dimension &VLENDataSet::extensionDim() const {
    return _extensionDim;
}

inline void VLENDataSet::flush() {
    if (valid() && H5Fflush(id(), H5F_SCOPE_LOCAL) < 0) {
        throw Exception("error when flushing HDF5 vlen data set with handle " + std::to_string(id()));
    }
}

inline void VLENDataSet::close() {
    auto pf = _parentFile.lock();
    if (pf) {
        if (!pf->closed() && valid() && H5Dclose(id()) < 0) {
            throw Exception("Error on closing HDF5 vlen data set");
        }
    }
}

inline VLENDataSet::~VLENDataSet() {
    try {
        close();
    } catch (const Exception &e) {
        std::cerr << "Unable to close hdf5 vlen data set: " << e.what() << std::endl;
    }
}

template<typename T>
inline void VLENDataSet::append(std::vector<std::vector<T>> &data) {
    if (!data.empty()) append({data.size()}, data.data());
}

template<typename T>
inline void VLENDataSet::append(const dimensions &dims, std::vector<T> *const data) {
    /*{
        std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("appending to vlen data set with data size = ({})", result.str());
    }*/
    {
        auto fs = getFileSpace();
        if (dims.size() != fs->ndim()) {
            throw Exception("Tried to append data with ndims=" + std::to_string(dims.size()) + " to set with ndims=" +
                            std::to_string(fs->ndim()));
        }
    }
    if (!_memorySpace) {
        _memorySpace = std::make_unique<DataSpace>(_parentFile, dims);
    } else {
        H5Sset_extent_simple(_memorySpace->id(), static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    dimensions currentExtent;
    {
        auto fs = getFileSpace();
        currentExtent = fs->dims();
    }
    dimensions offset;
    {
        offset.resize(dims.size());
        offset[_extensionDim] = currentExtent[_extensionDim];
    }
    dimensions newExtent(currentExtent.begin(), currentExtent.end());
    {
        newExtent[_extensionDim] += dims[_extensionDim];
        H5Dset_extent(id(), newExtent.data());
        {
            // todo log::trace("setting new extent to:");
            /*std::stringstream result;
            std::copy(newExtent.begin(), newExtent.end(), std::ostream_iterator<int>(result, ", "));
            log::trace("    size = {}", result.str());*/
        }
    }
    // todo log::trace("selecting hyperslab with:");
    {
        /*std::stringstream result;
        std::copy(offset.begin(), offset.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("    current extent = {}", result.str());*/
    }
    {
        /*std::stringstream result;
        std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>(result, ", "));
        log::trace("    size = {}", result.str());*/
    }
    std::vector<hvl_t> traj;
    {
        const auto n = dims[_extensionDim];
        // log::trace("appending n={} vlen entries", n);
        traj.reserve(n);
        for (auto i = 0; i < n; ++i) {
            auto val = &data[i];
            hvl_t entry{};
            entry.len = val->size();
            entry.p = val->data();
            traj.push_back(entry);
        }
    }
    {
        auto fileSpace = getFileSpace();
        H5Sselect_hyperslab(fileSpace->id(), H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
        H5Dwrite(id(), _memoryType.id(), _memorySpace->id(), fileSpace->id(), H5P_DEFAULT, traj.data());
    }
}

}