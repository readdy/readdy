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
 * @file traits.h
 * @brief << brief description >>
 * @author clonker
 * @date 04.09.17
 * @copyright BSD-3
 */

#pragma once

#include <string>
#include <array>
#include <vector>
#include <memory>

#include <H5Gpublic.h>
#include <H5Spublic.h>
#include <hdf5.h>

namespace h5rd {

class File;

class Filter;

class Group;

class DataSet;

class VLENDataSet;

template<typename T>
class STDDataSetType;

template<typename T>
class NativeDataSetType;

template<typename Container>
class Node;

class DataSpace;

class DataSetType;


using handle_id = hid_t;
using dimension = hsize_t;
using dimensions = std::vector<dimension>;
using group_info = H5G_info_t;
const static unsigned long long UNLIMITED_DIMS = H5S_UNLIMITED;

namespace util {

inline bool groupExists(hid_t hid, const std::string &name) {
    H5E_BEGIN_TRY
        {
            hid = H5Gopen(hid, name.c_str(), H5P_DEFAULT);
            if (hid > 0) {
                H5Gclose(hid);
            }
        }
    H5E_END_TRY
    return (hid > 0);
}

template<typename T>
struct n_dims {
    static constexpr std::size_t value = 0;
};

template<typename T>
struct n_dims<std::vector<T>> {
    static constexpr std::size_t value = 1 + n_dims<T>::value;
};

template<typename T>
struct n_dims<T *> {
    static constexpr std::size_t value = 1 + n_dims<T>::value;
};

template<typename T, std::size_t N>
struct n_dims<T[N]> {
    static constexpr std::size_t value = 1 + n_dims<T>::value;
};

template<typename T>
struct inner_type {
    using type = T;
};

template<typename T>
struct inner_type<std::vector<T>> {
    using type = typename inner_type<T>::type;
};

template<typename T>
struct inner_type<T *> {
    using type = typename inner_type<T>::type;
};

template<typename T, std::size_t N>
struct inner_type<T[N]> {
    using type = typename inner_type<T>::type;
};

template<typename T>
struct is_std_array : public std::false_type {
};

template<typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : public std::true_type {
};

template<typename... Ts>
struct make_void {
    typedef void type;
};
template<typename... Ts> using void_t = typename make_void<Ts...>::type;

}
}
