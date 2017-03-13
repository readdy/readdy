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
 * @file String_utils.h
 * @brief << brief description >>
 * @author clonker
 * @date 10.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include <hdf5.h>
#include <string>
#include <vector>
#include <algorithm>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

class WriteString {
public:
    WriteString(hid_t dataset, hid_t datatype, hid_t dataspace, hid_t memspace)
            : m_dataset(dataset), m_datatype(datatype), m_dataspace(dataspace), m_memspace(memspace), m_pos() {}

private:
    hid_t m_dataset;
    hid_t m_datatype;
    hid_t m_dataspace;
    hid_t m_memspace;
    hsize_t m_pos;

public:
    void operator()(std::vector<std::string>::value_type const &v) {
        // Select the file position, 1 record at position 'pos'
        hsize_t count[] = {1};
        hsize_t offset[] = {m_pos++};
        H5Sselect_hyperslab(m_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

        const char *s = v.c_str();
        H5Dwrite(m_dataset, m_datatype, m_memspace, m_dataspace, H5P_DEFAULT, &s);
    }
};

inline void writeVector(hid_t group, const std::string& dsName, std::vector<std::string> const &v) {
    hsize_t dims[] = {v.size()};
    hid_t dataspace = H5Screate_simple(sizeof(dims) / sizeof(*dims), dims, NULL);

    dims[0] = 1;
    hid_t memspace = H5Screate_simple(sizeof(dims) / sizeof(*dims), dims, NULL);

    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, H5T_VARIABLE);

    hid_t dataset = H5Dcreate1(group, dsName.c_str(), datatype, dataspace, H5P_DEFAULT);

    //
    // Select the "memory" to be written out - just 1 record.
    hsize_t offset[] = {0};
    hsize_t count[] = {1};
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    std::for_each(v.begin(), v.end(), WriteString(dataset, datatype, dataspace, memspace));

    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Tclose(datatype);
}


NAMESPACE_END(io)
NAMESPACE_END(readdy)
