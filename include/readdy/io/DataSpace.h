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
#include "H5Types.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

class DataSpaceHandle {
public:
    DataSpaceHandle(h5::handle_t hid) : hid(hid) {}

    ~DataSpaceHandle() {
        if (hid >= 0) {
            if (H5Sclose(hid) < 0) {
                log::error("failed to close data space {}!", hid);
                H5Eprint(H5Eget_current_stack(), stderr);
            }
        }
    }

    h5::handle_t operator*() const {
        return hid;
    }

private:
    h5::handle_t hid;
};

class DataSpace {
public:
    explicit DataSpace(h5::handle_t handle);

    explicit DataSpace(const std::vector<h5::dims_t> &dims, const std::vector<h5::dims_t> &maxDims = {});

    std::size_t ndim() const;

    std::vector<h5::dims_t> dims() const;

    h5::handle_t handle() const;

private:
    std::shared_ptr<DataSpaceHandle> _handle;
};

NAMESPACE_END(io)
NAMESPACE_END(readdy)

#include "bits/DataSpace_bits.h"
