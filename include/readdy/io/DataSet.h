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
 * @file DataSet.h
 * @brief << brief description >>
 * @author clonker
 * @date 04/01/2017
 * @copyright GNU Lesser General Public License v3.0
 */
#pragma once
#include "H5Types.h"
#include "Group.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

enum DataSetCompression {
    none = 0x0000, blosc = 0x0001
};

namespace blosc_compression {
void initialize();
void activate(hid_t plist, unsigned int* cd_values);
}

class DataSetHandle {
public:
    DataSetHandle(h5::handle_t handle = -1) : handle(handle) {}
    ~DataSetHandle() {
        if(handle >= 0) {
            if(H5Dclose(handle) < 0) {
                log::error("error on closing data set {}!", handle);
                H5Eprint(H5Eget_current_stack(), stderr);
            }
        }
    }

    h5::handle_t operator*() const {
        return handle;
    }
private:
    h5::handle_t handle;
};

template<typename T, bool VLEN=false, int compression=DataSetCompression::blosc>
class READDY_API DataSet {
public:

    DataSet(const std::string &name, const Group &group, const std::vector<h5::dims_t> &chunkSize,
            const std::vector<h5::dims_t> &maxDims);

    DataSet(const std::string &name, const Group &group, const std::vector<h5::dims_t> &chunkSize,
            const std::vector<h5::dims_t> &maxDims, DataSetType memoryType, DataSetType fileType);

    virtual ~DataSet();

    DataSet(DataSet&& rhs) = default;
    DataSet& operator=(DataSet&&) = delete;
    DataSet(const DataSet&) = delete;
    DataSet& operator=(const DataSet&) = delete;

    void close();

    template<bool no_vlen = !VLEN>
    void append(std::vector<T> &data, typename std::enable_if<no_vlen, bool>::type* = 0);

    template<bool no_vlen = !VLEN>
    void append(const std::vector<h5::dims_t> &dims, const T *const data, typename std::enable_if<no_vlen, bool>::type* = 0);

    template<bool vlen = VLEN>
    void append(std::vector<std::vector<T>> &data, typename std::enable_if<vlen, bool>::type* = 0);

    template<bool vlen = VLEN>
    void append(const std::vector<h5::dims_t> &dims, std::vector<T> *const data, typename std::enable_if<vlen, bool>::type* = 0);

    void flush();

private:
    const std::vector<h5::dims_t> maxDims;
    const Group group;
    std::shared_ptr<DataSetHandle> handle_ref;
    h5::dims_t extensionDim;
    h5::handle_t memorySpace;
    DataSetType memoryType {};
    DataSetType fileType {};
};

NAMESPACE_END(io)
NAMESPACE_END(readdy)

#include "bits/DataSet_bits.h"
