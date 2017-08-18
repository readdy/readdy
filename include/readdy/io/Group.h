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
 * @file Group.h
 * @brief << brief description >>
 * @author clonker
 * @date 04/01/2017
 * @copyright GNU Lesser General Public License v3.0
 */
#pragma once

#include <string>
#include <vector>
#include <readdy/common/macros.h>
#include "H5Types.h"
#include "DataSetType.h"
#include "Object.h"
#include "DataSet.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

class READDY_API GroupHandle : public ObjectHandle {
public:
    explicit GroupHandle(h5::h5_handle handle = -1) : ObjectHandle(handle) {}

    GroupHandle(const GroupHandle&) = default;
    GroupHandle& operator=(const GroupHandle&) = default;
    GroupHandle(GroupHandle&&) = default;
    GroupHandle& operator=(GroupHandle&&) = default;

    ~GroupHandle() {
        if (_handle >= 0) {
            close();
        }
    }

    void close() override {
        if (H5Gclose(_handle) < 0) {
            log::error("error on closing group!");
            H5Eprint(H5Eget_current_stack(), stderr);
        }
    }

};

class READDY_API Group : public Object {
    friend class File;

    friend class DataSet;

public:

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<T> &data) {
        write(dataSetName, {data.size()}, data.data());
    }

    void write(const std::string &dataSetName, const std::string &string);

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<h5::h5_dims> &dims, const T *data);

    Group createGroup(const std::string &path);

    std::vector<std::string> subgroups() const;

    std::vector<std::string> contained_data_sets() const;

    template<typename T>
    void read(const std::string &dataSetName, std::vector<T> &array);

    template<typename T>
    void read(const std::string &dataSetName, std::vector<T> &array, DataSetType memoryType, DataSetType fileType);

    Group subgroup(const std::string &name);

    h5::h5_group_info info() const;

    template<typename T>
    DataSet createDataSet(const std::string &name, const std::vector<h5::h5_dims> &chunkSize,
                          const std::vector<h5::h5_dims> &maxDims,
                          DataSetCompression compression = DataSetCompression::blosc) {
        return createDataSet(name, chunkSize, maxDims, STDDataSetType<T>(), NativeDataSetType<T>(), compression);
    }

    DataSet createDataSet(const std::string &name, const std::vector<h5::h5_dims> &chunkSize,
                          const std::vector<h5::h5_dims> &maxDims, const DataSetType &memoryType,
                          const DataSetType &fileType, DataSetCompression compression = DataSetCompression::blosc) {
        h5::h5_dims extensionDim;
        {
            std::stringstream result;
            std::copy(maxDims.begin(), maxDims.end(), std::ostream_iterator<int>(result, ", "));
            log::trace("creating data set with maxDims={}", result.str());
        }
        if (compression == DataSetCompression::blosc) {
            blosc_compression::initialize();
        }
        // validate and find extension dim
        {
            auto unlimited_it = std::find(maxDims.begin(), maxDims.end(), h5::UNLIMITED_DIMS);
            bool containsUnlimited = unlimited_it != maxDims.end();
            if (!containsUnlimited) {
                throw std::runtime_error("needs to contain unlimited_dims in some dimension to be extensible");
            }
            extensionDim = static_cast<h5::h5_dims>(std::distance(maxDims.begin(), unlimited_it));
            log::trace("found extension dim {}", extensionDim);
        }
        h5::h5_handle handle;
        {
            // set up empty data set
            std::vector<h5::h5_dims> dims(maxDims.begin(), maxDims.end());
            dims[extensionDim] = 0;
            DataSpace fileSpace(dims, maxDims);
            DataSetCreatePropertyList propertyList;
            propertyList.set_layout_chunked();
            propertyList.set_chunk(chunkSize);
            if (compression == DataSetCompression::blosc) propertyList.activate_blosc();

            auto _hid = H5Dcreate(hid(), name.c_str(), fileType.hid(),
                                  fileSpace.hid(), H5P_DEFAULT, propertyList.hid(), H5P_DEFAULT);
            if (_hid < 0) {
                log::error("Error on creating data set {}", _hid);
                H5Eprint(H5Eget_current_stack(), stderr);
                throw std::runtime_error("Error on creating data set " + std::to_string(_hid));
            } else {
                handle = _hid;
            }
        }
        DataSet ds(handle, memoryType, fileType);
        ds.extensionDim() = extensionDim;
        return ds;
    }

    template<typename T>
    VLENDataSet createVLENDataSet(const std::string &name, const std::vector<h5::h5_dims> &chunkSize,
                                  const std::vector<h5::h5_dims> &maxDims) {
        return createVLENDataSet(name, chunkSize, maxDims, STDDataSetType<T>(), NativeDataSetType<T>());
    }

    VLENDataSet createVLENDataSet(const std::string &name, const std::vector<h5::h5_dims> &chunkSize,
                                  const std::vector<h5::h5_dims> &maxDims, const DataSetType &memoryType,
                                  const DataSetType &fileType) {
        h5::h5_dims extensionDim;
        {
            std::stringstream result;
            std::copy(maxDims.begin(), maxDims.end(), std::ostream_iterator<int>(result, ", "));
            log::trace("creating vlen data set with maxDims={}", result.str());
        }
        VLENDataSetType vlenMemoryType(memoryType);
        VLENDataSetType vlenFileType(fileType);

        // validate and find extension dim
        {
            auto unlimited_it = std::find(maxDims.begin(), maxDims.end(), h5::UNLIMITED_DIMS);
            bool containsUnlimited = unlimited_it != maxDims.end();
            if (!containsUnlimited) {
                throw std::runtime_error("needs to contain unlimited_dims in some dimension to be extensible");
            }
            extensionDim = static_cast<h5::h5_dims>(std::distance(maxDims.begin(), unlimited_it));
            log::trace("found extension dim {}", extensionDim);
        }
        h5::h5_handle hid;
        {
            // set up empty data set
            std::vector<h5::h5_dims> dims(maxDims.begin(), maxDims.end());
            dims[extensionDim] = 0;
            DataSpace fileSpace(dims, maxDims);
            DataSetCreatePropertyList propertyList;
            propertyList.set_layout_chunked();
            propertyList.set_chunk(chunkSize);

            hid = H5Dcreate(Group::hid(), name.c_str(), vlenFileType.hid(),
                            fileSpace.hid(), H5P_DEFAULT, propertyList.hid(), H5P_DEFAULT);
            if (hid < 0) {
                log::error("Error on creating data set {}", hid);
                H5Eprint(H5Eget_current_stack(), stderr);
                throw std::runtime_error("Error on creating data set " + std::to_string(hid));
            }
        }
        VLENDataSet ds(hid, vlenMemoryType, vlenFileType);
        ds.extensionDim() = extensionDim;
        log::debug("set vlen extension dim to {}", ds.extensionDim());
        return ds;
    }

protected:

    Group();

    Group(h5::h5_handle handle, const std::string &);

    std::vector<std::string> sub_elements(H5O_type_t type) const;

    std::string path;
};

NAMESPACE_END(io)
NAMESPACE_END(readdy)

#include "bits/Group_bits.h"
