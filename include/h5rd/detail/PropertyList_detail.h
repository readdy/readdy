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
 * @file PropertyList_detail.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.09.17
 * @copyright BSD-3
 */

#pragma once

#include "../PropertyList.h"
#include "../Filter.h"

#include <iostream>


namespace h5rd {

inline PropertyList::~PropertyList() {
    try {
        close();
    } catch (const Exception &e) {
        std::cerr << "Unable to close hdf5 property list: " << e.what() << std::endl;
    }
}

inline void PropertyList::close() {
    if (valid()) {
        if (H5Pclose(id()) < 0) {
            throw Exception("Error on closing HDF5 property list");
        }
    }
}

inline PropertyList::PropertyList(handle_id cls_id, ParentFileRef parentFile) : super(std::move(parentFile)) {
    _hid = H5Pcreate(cls_id);
    if (!valid()) {
        std::stringstream ss;
        ss << "Failed to create property list! Reason: ";
        {
            auto pf = _parentFile.lock();
            if(!pf) {
                ss << "Parent file weak ptr was erased";
            } else {
                if(pf->closed()) {
                    ss << "Parent file was closed";
                } else {
                    if(_hid == H5I_INVALID_HID) {
                        ss << "HID was invalid id";
                    } else {
                        if (H5Iis_valid(_hid) <= 0) {
                            ss << "H5Iis_valid check failed.";
                        } else {
                            ss << "nothing wrong, should not happen...";
                        }
                    }
                }
            }
        }
        throw Exception(ss.str());
    }
}

inline LinkCreatePropertyList::LinkCreatePropertyList(ParentFileRef parentFile) : PropertyList(H5P_LINK_CREATE,
                                                                                               std::move(parentFile)) {}

inline void LinkCreatePropertyList::set_create_intermediate_group() {
    H5Pset_create_intermediate_group(id(), 1);
}

inline FileAccessPropertyList::FileAccessPropertyList(ParentFileRef parentFile) : PropertyList(H5P_FILE_ACCESS,
                                                                                               std::move(parentFile)) {}

inline void FileAccessPropertyList::set_close_degree_weak() {
    H5Pset_fclose_degree(id(), H5F_CLOSE_WEAK);
}

inline void FileAccessPropertyList::set_close_degree_semi() {
    H5Pset_fclose_degree(id(), H5F_CLOSE_SEMI);
}

inline void FileAccessPropertyList::set_close_degree_strong() {
    H5Pset_fclose_degree(id(), H5F_CLOSE_STRONG);
}

inline void FileAccessPropertyList::set_close_degree_default() {
    H5Pset_fclose_degree(id(), H5F_CLOSE_DEFAULT);
}

inline void FileAccessPropertyList::set_use_latest_libver() {
    H5Pset_libver_bounds(id(), H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
}

inline DataSetCreatePropertyList::DataSetCreatePropertyList(ParentFileRef parentFile)
                    : PropertyList(H5P_DATASET_CREATE, std::move(parentFile)) {}

inline void DataSetCreatePropertyList::set_layout_compact() {
    H5Pset_layout(id(), H5D_COMPACT);
}

inline void DataSetCreatePropertyList::set_layout_contiguous() {
    H5Pset_layout(id(), H5D_CONTIGUOUS);
}

inline void DataSetCreatePropertyList::set_layout_chunked() {
    H5Pset_layout(id(), H5D_CHUNKED);
}

inline void DataSetCreatePropertyList::set_chunk(const dimensions &chunk_dims) {
    H5Pset_chunk(id(), static_cast<int>(chunk_dims.size()), chunk_dims.data());
}

inline void DataSetCreatePropertyList::activate_filter(Filter *filter) {
    filter->activate(*this);
}


}