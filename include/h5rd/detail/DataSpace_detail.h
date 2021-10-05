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
 * @file DataSpace_detail.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include "../DataSpace.h"
#include <iostream>

inline std::size_t h5rd::DataSpace::ndim() const {
    const auto n = H5Sget_simple_extent_ndims(id());
    if (n < 0) {
        throw Exception("Failed to retrieve ndims for data space");
    }
    return static_cast<std::size_t>(n);
}

inline h5rd::dimensions h5rd::DataSpace::dims() const {
    if (!valid()) {
        throw Exception("Tried requesting dims from invalid data space");
    }
    dimensions result;
    result.resize(ndim());
    if (H5Sget_simple_extent_dims(id(), result.data(), nullptr) < 0) {
        throw Exception("Failed to retrieve dims for data space");
    }
    return result;
}

inline h5rd::dimensions h5rd::DataSpace::maxDims() const {
    if (!valid()) {
        throw Exception("Tried requesting maxDims from invalid data space");
    }
    dimensions result;
    result.resize(ndim());
    if (H5Sget_simple_extent_dims(id(), nullptr, result.data()) < 0) {
        throw Exception("Failed to retrieve max dims for data space");
    }
    return result;
}

inline void h5rd::DataSpace::close() {
    auto pf = _parentFile.lock();
    if (pf) {
        if (!pf->closed() && valid()) {
            if (H5Sclose(id()) < 0) {
                throw Exception("Error on closing HDF5 data space");
            }
        }
    }
}

inline h5rd::DataSpace::~DataSpace() {
    try {
        close();
    } catch (const Exception &e) {
        std::cerr << "Unable to close hdf5 data space: " << e.what() << std::endl;
    }
}

inline h5rd::DataSpace::DataSpace(ParentFileRef parentFile, const dimensions &dims, const dimensions &maxDims)
        : SubObject(std::move(parentFile)) {
    if (maxDims.empty()) {
        _hid = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
    } else {
        _hid = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), maxDims.data());
    }
    if (_hid < 0) {
        throw Exception("Error on creating data space!");
    }
}

// inline h5rd::DataSpace::DataSpace() : Object(nullptr) {}

inline h5rd::DataSpace::DataSpace(ParentFileRef parentFile, h5rd::handle_id handle) : SubObject(std::move(parentFile)) {
    _hid = handle;
}
