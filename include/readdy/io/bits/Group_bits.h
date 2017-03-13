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
#include <readdy/io/DataSetType.h>
#include <readdy/io/File.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

inline Group::Group(h5::handle_t handle, const std::string &path)
        : handle(handle), path(path) {}

inline Group::Group() : Group(-1, "/") {}

inline void Group::write(const std::string &dataSetName, const std::string &string) {
    auto stdstr = STDDataSetType<std::string>();
    auto nativestr = NativeDataSetType<std::string>();
    H5Tset_cset(stdstr.tid->tid, H5T_CSET_UTF8);
    H5Tset_cset(nativestr.tid->tid, H5T_CSET_UTF8);
    H5LTmake_dataset_string(handle, dataSetName.c_str(), string.c_str());
}

inline Group Group::createGroup(const std::string &path) {
    if(util::groupExists(*this, path)) {
        return Group(H5Gopen(this->handle, path.c_str(), H5P_DEFAULT), path);
    } else {
        auto plist = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(plist, 1);
        auto handle = H5Gcreate(this->handle, path.c_str(), plist, H5P_DEFAULT, H5P_DEFAULT);
        return Group(handle, path);
    }
}

inline h5::handle_t Group::getHandle() const {
    return handle;
}

template<>
inline void Group::write<std::string>(const std::string &dataSetName, const std::vector<std::string> &data) {
    writeVector(handle, dataSetName, data);
}

template<>
inline void Group::write<short>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const short *data) {
    H5LTmake_dataset_short(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void Group::write<int>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const int *data) {
    H5LTmake_dataset_int(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void Group::write<long>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const long *data) {
    H5LTmake_dataset_long(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void Group::write<float>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const float *data) {
    H5LTmake_dataset_float(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void Group::write<double>(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const double *data) {
    H5LTmake_dataset_double(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

NAMESPACE_END(io)
NAMESPACE_END(readdy)
