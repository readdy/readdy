/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file String_utils.h
 * @brief << brief description >>
 * @author clonker
 * @date 10.03.17
 * @copyright BSD-3
 */

#pragma once

#include <string>
#include <vector>
#include <algorithm>

#include <H5Ipublic.h>
#include <H5Spublic.h>
#include <H5Dpublic.h>
#include <H5Tpublic.h>
#include <H5Ppublic.h>

namespace h5rd {
namespace util {
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
        H5Sselect_hyperslab(m_dataspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);

        const char *s = v.c_str();
        H5Dwrite(m_dataset, m_datatype, m_memspace, m_dataspace, H5P_DEFAULT, &s);
    }
};

template<typename T>
inline void writeVector(hid_t group, const std::string &dsName, std::vector<T> const &v) {};

template<>
inline void writeVector(hid_t group, const std::string &dsName, std::vector<std::string> const &v) {
    hsize_t dims[] = {v.size()};
    hid_t dataspace = H5Screate_simple(sizeof(dims) / sizeof(*dims), dims, nullptr);

    dims[0] = 1;
    hid_t memspace = H5Screate_simple(sizeof(dims) / sizeof(*dims), dims, nullptr);

    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, H5T_VARIABLE);

    hid_t dataset = H5Dcreate1(group, dsName.c_str(), datatype, dataspace, H5P_DEFAULT);

    //
    // Select the "memory" to be written out - just 1 record.
    hsize_t offset[] = {0};
    hsize_t count[] = {1};
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    std::for_each(v.begin(), v.end(), WriteString(dataset, datatype, dataspace, memspace));

    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Tclose(datatype);
}
}
}
