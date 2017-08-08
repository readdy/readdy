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
 * @file PropertyList.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.05.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include "Object.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)


namespace blosc_compression {
void initialize();
void activate(hid_t plist, unsigned int* cd_values);
}


class PropertyListHandle : public ObjectHandle {
public:
    explicit PropertyListHandle(h5::handle_t handle) : ObjectHandle(handle) {}

    PropertyListHandle(const PropertyListHandle&) = default;
    PropertyListHandle& operator=(const PropertyListHandle&) = default;
    PropertyListHandle(PropertyListHandle&&) = default;
    PropertyListHandle& operator=(PropertyListHandle&&) = default;

    ~PropertyListHandle() override {
        if(_handle >= 0) close();
    }

    void close() override {
        if(H5Pclose(_handle) < 0) {
            log::error("error on closing property list");
            H5Eprint(H5Eget_current_stack(), stderr);
        }
    }
};

class PropertyList : public Object {

protected:
    explicit PropertyList(h5::handle_t cls_id) : Object(std::make_shared<PropertyListHandle>(cls_id)) {}

};

class DataSetCreatePropertyList : public PropertyList {
public:
    DataSetCreatePropertyList() : PropertyList(H5Pcreate(H5P_DATASET_CREATE)) {}

    void set_layout_compact() {
        H5Pset_layout(hid(), H5D_COMPACT);
    }

    void set_layout_contiguous() {
        H5Pset_layout(hid(), H5D_CONTIGUOUS);
    }

    void set_layout_chunked() {
        H5Pset_layout(hid(), H5D_CHUNKED);
    }

    void set_chunk(const std::vector<h5::dims_t> &chunk_dims) {
        H5Pset_chunk(hid(), static_cast<int>(chunk_dims.size()), chunk_dims.data());
    }

    void activate_blosc() {
        unsigned int cd_values[7];
        blosc_compression::activate(hid(), cd_values);
    }


};

NAMESPACE_END(io)
NAMESPACE_END(readdy)
