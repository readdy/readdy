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

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

class READDY_API GroupHandle {
public:
    GroupHandle(h5::handle_t handle = -1) : handle(handle) {}

    ~GroupHandle() {
        if (handle >= 0) {
            if (H5Gclose(handle) < 0) {
                log::error("error when closing group {}", handle);
                H5Eprint(H5Eget_current_stack(), stderr);
            }
        }
    }

    void set(h5::handle_t handle) {
        GroupHandle::handle = handle;
    }

    h5::handle_t operator*() {
        return handle;
    }

private:
    h5::handle_t handle{-1};
};

class READDY_API Group {
    friend class File;

    template<typename T, bool VLEN, int compression>
    friend
    class DataSet;

public:
    using handle_ref = std::shared_ptr<GroupHandle>;

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<T> &data) {
        write(dataSetName, {data.size()}, data.data());
    }

    void write(const std::string &dataSetName, const std::string &string);

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const T *data);

    Group createGroup(const std::string &path);

    h5::handle_t getHandle() const;

    std::vector<std::string> subgroups() const;

    Group subgroup(const std::string &name);

    h5::group_info_t info() const;

protected:

    Group();

    Group(h5::handle_t handle, const std::string &);

    handle_ref handle;
    std::string path;
};

NAMESPACE_END(io)
NAMESPACE_END(readdy)

#include "bits/Group_bits.h"
