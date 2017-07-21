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
 * @file Group_bits.h
 * @brief << brief description >>
 * @author clonker
 * @date 04/01/2017
 * @copyright GNU Lesser General Public License v3.0
 */
#pragma once
#include "../Group.h"
#include "Util_bits.h"
#include "String_utils.h"
#include <hdf5_hl.h>
#include "../DataSetType.h"
#include "../File.h"
#include <H5Lpublic.h>
#include <readdy/io/DataSetType.h>
#include "../DataSpace.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

inline Group::Group(h5::handle_t handle, const std::string &path)
        : Object(std::make_shared<GroupHandle>(handle)), path(path) {}

inline Group::Group() : Group(-1, "/") {}

inline void Group::write(const std::string &dataSetName, const std::string &string) {
    auto stdstr = STDDataSetType<std::string>();
    auto nativestr = NativeDataSetType<std::string>();
    H5Tset_cset(stdstr.hid(), H5T_CSET_UTF8);
    H5Tset_cset(nativestr.hid(), H5T_CSET_UTF8);
    if(H5LTmake_dataset_string(**handle, dataSetName.c_str(), string.c_str()) < 0) {
        log::warn("there was a problem with writing {} into a hdf5 file.", string);
    }
}

inline Group Group::createGroup(const std::string &path) {
    if(util::groupExists(*this, path)) {
        return Group(H5Gopen(**(this->handle), path.c_str(), H5P_DEFAULT), path);
    } else {
        auto plist = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(plist, 1);
        auto handle = H5Gcreate(**(this->handle), path.c_str(), plist, H5P_DEFAULT, H5P_DEFAULT);
        return Group(handle, path);
    }
}

inline h5::group_info_t Group::info() const {
    h5::group_info_t info;
    if(H5Gget_info(hid(), &info) < 0) {
        log::critical("Failure to get group info!");
        H5Eprint(H5Eget_current_stack(), stderr);
    }
    return info;
}

inline std::vector<std::string> Group::sub_elements(H5O_type_t type) const {
    std::vector<std::string> result;
    auto group_info = info();
    result.reserve(group_info.nlinks);
    for(std::size_t i=0; i < group_info.nlinks; ++i) {
        H5O_info_t oinfo;
        H5Oget_info_by_idx(hid(), ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo, H5P_DEFAULT);
        if(oinfo.type == type) {
            auto size = 1+ H5Lget_name_by_idx (hid(), ".", H5_INDEX_NAME, H5_ITER_INC, i, NULL, 0, H5P_DEFAULT);
            if(size < 0) {
                H5Eprint(H5Eget_current_stack(), stderr);
            }

            std::vector<char> c_string (static_cast<std::size_t>(size));
            H5Lget_name_by_idx (hid(), ".", H5_INDEX_NAME, H5_ITER_INC, i, c_string.data(), (std::size_t) size, H5P_DEFAULT);
            std::string label (c_string.data());

            result.push_back(std::move(label));
        }
    }
    return result;
}

inline std::vector<std::string> Group::subgroups() const {
    return sub_elements(H5O_TYPE_GROUP);
}

inline std::vector<std::string> Group::contained_data_sets() const {
    return sub_elements(H5O_TYPE_DATASET);
}

inline Group Group::subgroup(const std::string& name) {
    auto gid = H5Gopen(hid(), name.data(), H5P_DEFAULT);
    if(gid < 0) {
        log::critical("could not open subgroup {}", name);
        H5Eprint(H5Eget_current_stack(), stderr);
    }
    return Group(gid, name);
}

template<typename T>
inline void Group::read(const std::string& name, std::vector<T> &array) {
    read(name, array, STDDataSetType<T>(), NativeDataSetType<T>());
}

template<typename T>
inline void Group::read(const std::string& ds_name, std::vector<T> &array, DataSetType memoryType, DataSetType fileType) {
    blosc_compression::initialize();

    const auto n_array_dims = 1 + util::n_dims<T>::value;
    auto hid = H5Dopen2(**handle, ds_name.data(), H5P_DEFAULT);

    DataSpace memorySpace (H5Dget_space(hid));

    const auto ndim = memorySpace.ndim();

    //if(ndim != n_array_dims) {
    //    log::error("wrong dimensionality: {} != {}", ndim, n_array_dims);
    //    throw std::invalid_argument("wrong dimensionality when attempting to read a data set");
    //}

    const auto dims = memorySpace.dims();
    std::size_t required_length = 1;
    for(const auto dim : dims) {
        log::trace("dim len = {}", dim);
        required_length *= dim;
    }
    log::trace("required length = {}", required_length);
    array.resize(required_length);

    auto result = H5Dread(hid, memoryType.hid(), H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data());

    if(result < 0) {
        log::trace("Failed reading result!");
        H5Eprint(H5Eget_current_stack(), stderr);
    }

    H5Dclose(hid);

    //for(std::size_t d = 0; d < ndim-1; ++d) {
    //    for(auto& sub_arr : array) {
    //        sub_arr.resize(dims[1]);
    //    }
    //}

    // todo reshape array to dims
}

template<>
inline void Group::write<std::string>(const std::string &dataSetName, const std::vector<std::string> &data) {
    writeVector(hid(), dataSetName, data);
}

template<>
inline void Group::write<short>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const short *data) {
    H5LTmake_dataset_short(hid(), dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void Group::write<int>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const int *data) {
    H5LTmake_dataset_int(hid(), dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void Group::write<long>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const long *data) {
    H5LTmake_dataset_long(hid(), dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void Group::write<float>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const float *data) {
    H5LTmake_dataset_float(hid(), dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void Group::write<double>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const double *data) {
    H5LTmake_dataset_double(hid(), dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

NAMESPACE_END(io)
NAMESPACE_END(readdy)
