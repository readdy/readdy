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
#ifndef READDY_MAIN_DATASET_H
#define READDY_MAIN_DATASET_H

#include "H5Types.h"
#include "Group.h"

namespace readdy {
namespace io {

template<typename T, bool VLEN=false>
class READDY_API DataSet {
public:

    DataSet(const std::string &name, const Group &group, const std::vector<h5::dims_t> &chunkSize,
            const std::vector<h5::dims_t> &maxDims);

    DataSet(const std::string &name, const Group &group, const std::vector<h5::dims_t> &chunkSize,
            const std::vector<h5::dims_t> &maxDims, DataSetType memoryType, DataSetType fileType);

    virtual ~DataSet();

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
    h5::dims_t extensionDim;
    h5::handle_t dataSetHandle;
    h5::handle_t memorySpace;
    DataSetType memoryType {};
    DataSetType fileType {};
};

}
}

#include "bits/DataSet_bits.h"

#endif //READDY_MAIN_DATASET_H
