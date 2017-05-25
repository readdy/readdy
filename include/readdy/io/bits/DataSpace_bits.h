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
 * @file DataSpace_bits.h
 * @brief << brief description >>
 * @author clonker
 * @date 25.05.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "../DataSpace.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

inline DataSpace::DataSpace(h5::handle_t handle) {
    _handle = std::make_shared<DataSpaceHandle>(handle);
}

inline DataSpace::DataSpace(const std::vector<h5::dims_t> &dims, const std::vector<h5::dims_t> &maxDims) {
    h5::handle_t hid;
    if (maxDims.empty()) {
        hid = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
    } else {
        hid = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), maxDims.data());
    }
    if (hid < 0) {
        log::error("error on creating data space!");
        H5Eprint(H5Eget_current_stack(), stderr);
    } else {
        _handle = std::make_shared<DataSpaceHandle>(hid);
    }
}

inline std::size_t DataSpace::ndim() const {
    const auto n = H5Sget_simple_extent_ndims(**_handle);
    if (n < 0) {
        log::error("failed to retrieve ndims");
        H5Eprint(H5Eget_current_stack(), stderr);
    }
    return static_cast<std::size_t>(n);
}

inline h5::handle_t DataSpace::handle() const {
    return **_handle;
}

inline std::vector<h5::dims_t> DataSpace::dims() const {
    std::vector<hsize_t> result;
    result.resize(ndim());
    if (H5Sget_simple_extent_dims(**_handle, result.data(), NULL) < 0) {
        log::error("failed to get dims!");
        H5Eprint(H5Eget_current_stack(), stderr);
    }
    return result;
}

NAMESPACE_END(io)
NAMESPACE_END(readdy)
