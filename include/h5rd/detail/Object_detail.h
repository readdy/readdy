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
 * @file Object_detail.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include "../Object.h"

inline const h5rd::Object::ParentFileRef &h5rd::Object::parentFile() const {
    return _parentFile;
}

inline bool h5rd::Object::closed() const {
    return _closed;
}

inline bool h5rd::Object::valid() const {
    auto pf = _parentFile.lock();
    if(pf) {
        return !pf->closed() && _hid != H5I_INVALID_HID && H5Iis_valid(_hid) > 0;
    }
    return false;
}

inline h5rd::handle_id h5rd::Object::id() const {
    return _hid;
}

inline h5rd::Object::Object(h5rd::Object &&rhs) noexcept
        : _hid(std::move(rhs._hid)), _parentFile(std::move(rhs._parentFile)), _closed(std::move(rhs._closed)) {
    rhs._closed = true;
}

inline h5rd::Object &h5rd::Object::operator=(h5rd::Object &&rhs) noexcept {
    _hid = rhs._hid;
    _parentFile = std::move(rhs._parentFile);
    _closed = rhs._closed;
    rhs._closed = true;
    return *this;
}

inline h5rd::SubObject::SubObject(ParentFileRef parentFile) : Object() {
    _parentFile = std::move(parentFile);
}
