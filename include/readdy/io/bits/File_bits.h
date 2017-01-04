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
 * @file File_bits.h
 * @brief << brief description >>
 * @author clonker
 * @date 04/01/2017
 * @copyright GNU Lesser General Public License v3.0
 */
#ifndef READDY_MAIN_FILE_BITS_H
#define READDY_MAIN_FILE_BITS_H

#include "../File.h"

namespace readdy {
namespace io {


inline unsigned getFlagValue(const File::Flag &flag) {
    switch (flag) {
        case File::Flag::READ_ONLY:
            return H5F_ACC_RDONLY;
        case File::Flag::READ_WRITE:
            return H5F_ACC_RDWR;
        case File::Flag::OVERWRITE:
            return H5F_ACC_TRUNC;
        case File::Flag::FAIL_IF_EXISTS:
            return H5F_ACC_EXCL;
        case File::Flag::CREATE_NON_EXISTING:
            return H5F_ACC_CREAT;
        case File::Flag::DEFAULT:
            return H5F_ACC_RDWR | H5F_ACC_CREAT | H5F_ACC_TRUNC;
    }
}

inline void File::close() {
    flush();
    if (root.handle >= 0 && H5Fclose(root.handle) < 0) {
        throw std::runtime_error("error when closing HDF5 file \"" + path_ + "\" with handle "
                                 + std::to_string(root.handle));
    } else {
        // we are now in closed state, set the handle to -1
        root.handle = -1;
    }
}

inline void File::flush() {
    if (root.handle >= 0 && H5Fflush(root.handle, H5F_SCOPE_LOCAL) < 0) {
        throw std::runtime_error("error when flushing HDF5 file \"" + path_ + "\" with handle "
                                 + std::to_string(root.handle));
    }
}

inline Group File::createGroup(const std::string &path) {
    auto handle = H5Gcreate(root.handle, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    return Group(handle, path);
}

inline File::~File() {
    close();
}

inline File::File(const std::string &path, const Action &action, const Flag &flag) : File(path, action,
                                                                                   std::vector<Flag>{flag}) {}

inline File::File(const std::string &path, const File::Action &action, const std::vector<File::Flag> &flags)
        : path_(path), root() {
    unsigned flag = 0x0000u;
    for (const auto &f : flags) {
        flag = flag | getFlagValue(f);
    }
    switch (action) {
        case Action::CREATE: {
            root.handle = H5Fcreate(path.c_str(), flag, H5P_DEFAULT, H5P_DEFAULT);
            break;
        }
        case Action::OPEN: {
            root.handle = H5Fopen(path.c_str(), flag, H5P_DEFAULT);
            break;
        }
    }
}

inline void File::write(const std::string &dataSetName, const std::string &data) {
    root.write(dataSetName, data);
}

inline const Group &File::getRootGroup() const {
    return root;
}


}
}
#endif //READDY_MAIN_FILE_BITS_H
