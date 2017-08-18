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
 * @file DataSpace.h
 * @brief << brief description >>
 * @author clonker
 * @date 25.05.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/common.h>
#include "Object.h"
#include "H5Types.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

class DataSpaceHandle : public ObjectHandle {
public:
    explicit DataSpaceHandle(h5::h5_handle hid) : ObjectHandle(hid) {}
    
    DataSpaceHandle(const DataSpaceHandle&) = default;
    DataSpaceHandle& operator=(const DataSpaceHandle&) = default;
    DataSpaceHandle(DataSpaceHandle&&) = default;
    DataSpaceHandle& operator=(DataSpaceHandle&&) = default;
    
    ~DataSpaceHandle() {
        if (_handle >= 0) {
            close();
        }
    }

    void close() override {
        if (H5Sclose(_handle) < 0) {
            log::error("failed to close data space {}!", _handle);
            H5Eprint(H5Eget_current_stack(), stderr);
        }
    }
};

class DataSpace : public Object {
public:
    explicit DataSpace(h5::h5_handle handle);

    explicit DataSpace(const std::vector<h5::h5_dims> &dims, const std::vector<h5::h5_dims> &maxDims = {});

    std::size_t ndim() const;

    std::vector<h5::h5_dims> dims() const;

    std::vector<h5::h5_dims> max_dims() const;
};

NAMESPACE_END(io)
NAMESPACE_END(readdy)

#include "bits/DataSpace_bits.h"
