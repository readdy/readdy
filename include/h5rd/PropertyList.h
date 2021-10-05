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
 * @file PropretySet.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include "Object.h"

namespace h5rd {

class PropertyList : public SubObject {
    using super = SubObject;
public:
    explicit PropertyList(handle_id cls_id, ParentFileRef parentFile);

    PropertyList(const PropertyList &) = delete;

    PropertyList(PropertyList &&) noexcept = delete;

    PropertyList &operator=(const PropertyList &) = delete;

    PropertyList &operator=(PropertyList &&) noexcept = delete;

    ~PropertyList() override;

    void close() override;
};

class LinkCreatePropertyList : public PropertyList {
public:
    explicit LinkCreatePropertyList(ParentFileRef parentFile);

    void set_create_intermediate_group();
};

class FileAccessPropertyList : public PropertyList {
public:
    explicit FileAccessPropertyList(ParentFileRef parentFile);

    void set_close_degree_weak();

    void set_close_degree_semi();

    void set_close_degree_strong();

    void set_close_degree_default();

    void set_use_latest_libver();
};

class DataSetCreatePropertyList : public PropertyList {
public:
    explicit DataSetCreatePropertyList(ParentFileRef parentFile);

    void set_layout_compact();

    void set_layout_contiguous();

    void set_layout_chunked();

    void set_chunk(const dimensions &chunk_dims);

    void activate_filter(Filter *filter);

};

}

#include "detail/PropertyList_detail.h"
