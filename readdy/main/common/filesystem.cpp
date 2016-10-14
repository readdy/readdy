/**
 * << detailed description >>
 *
 * @file filesystem.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 14.10.16
 */

#include <stdio.h>
#include <sys/stat.h>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#include <system_error>
#define GetCurrentDir getcwd
#endif

#include <readdy/common/filesystem.h>
#include <readdy/common/make_unique.h>
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
    struct stat s;
    if (stat(path.c_str(), &s) == 0) {
        if (s.st_mode & S_IFREG) {
            return true;
        }
    }
    return false;
}

bool is_directory(const std::string &path) {
    struct stat s;
    if (stat(path.c_str(), &s) == 0) {
        if (s.st_mode & S_IFDIR) {
            return true;
        }
    }
    return false;
}

struct dir::Impl {
    tinydir_dir dir;
};


dir::dir(const std::string &path) : pimpl(std::make_unique<Impl>()) {
    if(exists(path)) {
        if(is_directory(path)) {
            if (tinydir_open(&(pimpl->dir), path.c_str()) == -1) {
                throw std::runtime_error("Error getting file " + path);
            }
        } else {
            std::runtime_error("file " + path + " did exist but was no directory");
        }
    } else {
        throw std::runtime_error("file " + path + " did not exist");
    }
}

dir::~dir() {
    tinydir_close(&(pimpl->dir));
}

bool dir::has_next() const {
    return static_cast<bool>(pimpl->dir.has_next);
}

std::string dir::next() {
    tinydir_file file;
    if(tinydir_readfile(&(pimpl->dir), &file) == -1) {
        throw std::runtime_error("error on getting file");
    }
    return std::string(file.name);
}
}
}
}