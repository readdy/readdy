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