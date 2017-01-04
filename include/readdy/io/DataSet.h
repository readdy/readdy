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

template<typename T>
class DataSet {
public:

    DataSet(const std::string &name, const Group &group, const std::vector<h5::dims_t> &chunkSize,
            const std::vector<h5::dims_t> &maxDims);

    ~DataSet();

    void close();

    template<typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void append(const std::vector<T> &data) {
        append({1, data.size()}, data.data());
    }

    void append(const std::vector<h5::dims_t> &dims, const T *data);

private:
    const std::vector<h5::dims_t> maxDims;
    const Group group;
    h5::dims_t extensionDim;
    h5::handle_t handle;
    h5::handle_t memorySpace;
    STDDataSetType<T> stdType {};
};

}
}

#include "bits/DataSet_bits.h"

#endif //READDY_MAIN_DATASET_H
