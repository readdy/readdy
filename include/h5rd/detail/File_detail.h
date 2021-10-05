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
 * @file File_detail.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include <iostream>

#include "../File.h"

namespace {

inline int convertFlag(const h5rd::File::Flag &flag) {
    switch (flag) {
        case h5rd::File::Flag::READ_ONLY:
            return H5F_ACC_RDONLY;
        case h5rd::File::Flag::READ_WRITE:
            return H5F_ACC_RDWR;
        case h5rd::File::Flag::OVERWRITE:
            return H5F_ACC_TRUNC;
        case h5rd::File::Flag::FAIL_IF_EXISTS:
            return H5F_ACC_EXCL;
        case h5rd::File::Flag::CREATE_NON_EXISTING:
            return H5F_ACC_CREAT;
        case h5rd::File::Flag::DEFAULT:
            return H5F_ACC_RDWR | H5F_ACC_CREAT | H5F_ACC_TRUNC;
    }
    throw std::logic_error("Unknown flag in convertFlag()!");
}
}

inline h5rd::File::File(const std::string &path, const Action &action, const Flag &flag)
        : File(path, action, Flags{flag}) {}

inline h5rd::File::File(const std::string &path, const Action &action, const Flags &flags)
        : Object(), path(path), action(action), flags(flags) {}

inline void h5rd::File::flush() {
    if (H5Fflush(_hid, H5F_SCOPE_LOCAL) < 0) {
        throw Exception("error when flushing HDF5 file \"" + path + "\"");
    }
}

inline h5rd::File::~File() {
    try {
        close();
    } catch (const Exception &e) {
        std::cerr << "Unable to close hdf5 file: " << e.what() << std::endl;
    }
}

inline void h5rd::File::close() {
    if (!closed() && _hid != H5I_INVALID_HID && H5Iis_valid(_hid) > 0) {
        if (H5Fclose(id()) < 0) {
            throw Exception("Error on closing HDF5 file \"" + path + "\"");
        }
        _closed = true;
    }
}

inline std::shared_ptr<h5rd::File> h5rd::File::open(const std::string &path, const h5rd::File::Flag &flag) {
    auto f = std::shared_ptr<h5rd::File>(new File(path, h5rd::File::Action::OPEN, flag));
    f->_parentFile = f->getptr();
    setUp(f);
    return f;
}

inline std::shared_ptr<h5rd::File> h5rd::File::open(const std::string &path, const h5rd::File::Flags &flags) {
    auto f = std::shared_ptr<h5rd::File>(new File(path, h5rd::File::Action::OPEN, flags));
    f->_parentFile = f->getptr();
    setUp(f);
    return f;
}

inline std::shared_ptr<h5rd::File> h5rd::File::create(const std::string &path, const h5rd::File::Flag &flag) {
    auto f = std::shared_ptr<h5rd::File>(new File(path, h5rd::File::Action::CREATE, flag));
    f->_parentFile = f->getptr();
    setUp(f);
    return f;
}

inline std::shared_ptr<h5rd::File> h5rd::File::create(const std::string &path, const h5rd::File::Flags &flags) {
    auto f = std::shared_ptr<h5rd::File>(new File(path, h5rd::File::Action::CREATE, flags));
    f->_parentFile = f->getptr();
    setUp(f);
    return f;
}

inline void h5rd::File::setUp(std::shared_ptr<File> file) {
    unsigned flag = 0x0000u;
    for (const auto &f : file->flags) {
        flag = flag | convertFlag(f);
    }
    FileAccessPropertyList fapl(file);
    fapl.set_close_degree_strong();
    fapl.set_use_latest_libver();
    handle_id val = 0;
    switch (file->action) {
        case Action::CREATE: {
            val = H5Fcreate(file->path.c_str(), flag, H5P_DEFAULT, fapl.id());
            break;
        }
        case Action::OPEN: {
            val = H5Fopen(file->path.c_str(), flag, fapl.id());
            break;
        }
    }
    if (val < 0) {
        throw Exception("Failed on opening/creating file " + file->path);
    }
    file->_hid = val;
}

inline std::shared_ptr<h5rd::Object> h5rd::File::getptr() {
    return shared_from_this();
}

inline h5rd::Object::ParentFileRef h5rd::File::ref() const {
    return parentFile();
}
