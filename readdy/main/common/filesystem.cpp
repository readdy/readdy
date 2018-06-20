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
 * @file filesystem.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 14.10.16
 */

#include <cstdio>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#include <system_error>
#include <readdy/common/logging.h>

#define GetCurrentDir getcwd
#endif

#include "readdy/common/filesystem.h"
#include "readdy/common/tinydir.h"

namespace readdy {
namespace util {
namespace fs {

std::string current_path()  {
    char cCurrentPath[FILENAME_MAX];

    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
    {
        throw std::system_error(errno, std::system_category());
    }
    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0';
    return std::string(cCurrentPath);
}

bool exists(const std::string& path) {
    return access( path.c_str(), F_OK ) != -1;
}

bool is_file(const std::string &path) {
    return !is_directory(path);
}

bool is_directory(const std::string &path) {
    tinydir_file file;
    if(tinydir_file_open(&file, path.c_str()) == -1) {
        log::error("error on opening {}", path);
        throw std::runtime_error("error on opening " + path);
    }
    return file.is_dir != 0;
}

struct dir_iterator::Impl {
    std::unique_ptr<tinydir_dir> dir = std::make_unique<tinydir_dir>();
    std::size_t n = 0;
};


dir_iterator::dir_iterator(const std::string &path) : pimpl(std::make_unique<Impl>()) {
    log::debug("getting dir iterator for path=\"{}\"", path);
    if(exists(path)) {
        if(is_directory(path)) {
            if (tinydir_open_sorted(pimpl->dir.get(), path.c_str()) == -1) {
                throw std::runtime_error("Error getting file " + path);
            }
        } else {
            std::runtime_error("file " + path + " did exist but was no directory");
        }
    } else {
        throw std::runtime_error("file " + path + " did not exist");
    }
}

dir_iterator::~dir_iterator() {
    tinydir_close(pimpl->dir.get());
}

bool dir_iterator::has_next() const {
    return pimpl->n < pimpl->dir->n_files;
}

std::string dir_iterator::next() {
    tinydir_file file;
    if(tinydir_readfile_n(pimpl->dir.get(), &file, pimpl->n++) == -1) {
        throw std::runtime_error("error on getting file");
    }
    return std::string(file.path);
}

dir_iterator::dir_iterator(dir_iterator && rhs) noexcept = default;

dir_iterator &dir_iterator::operator=(dir_iterator && rhs) noexcept = default;

std::string dir_iterator::base_name() const {
    const std::string path{pimpl->dir->path};
    auto lastSlash = path.find_last_of("/\\");
    if (lastSlash == std::string::npos) {
        return path;
    }
    return path.substr(lastSlash + 1);
}

}
}
}