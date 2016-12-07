/**
 * << detailed description >>
 *
 * @file File.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 31.08.16
 */

#include "readdy/io/File.h"

#include <stdexcept>

#include <hdf5.h>
#include <hdf5_hl.h>

namespace readdy {
namespace io {

unsigned getFlagValue(const File::Flag &flag) {
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
    }
}

Group::Group(Group::handle_t handle) : handle(handle){

}

void File::close() {
    if (self.handle >= 0) {
        flush();
        H5Fclose(self.handle);
    }
}

void File::flush() {
    if (self.handle >= 0 && H5Fflush(self.handle, H5F_SCOPE_LOCAL) < 0) {
        throw std::runtime_error("flushing HDF5 file: " + path_);
    }
}

File::~File() {
    close();
}

File::File(const std::string &path, const Action &action, const Flag &flag) : File(path, action,
                                                                                   std::vector<Flag>{flag}) {}

File::File(const std::string &path, const File::Action &action, const std::vector<File::Flag> &flags)
        : path_(path), self({}) {
    unsigned flag = 0x0000u;
    for (const auto &f : flags) {
        flag = flag | getFlagValue(f);
    }
    switch (action) {
        case Action::CREATE: {
            self.handle = H5Fcreate(path.c_str(), flag, H5P_DEFAULT, H5P_DEFAULT);
            break;
        }
        case Action::OPEN: {
            self.handle = H5Fopen(path.c_str(), flag, H5P_DEFAULT);
            break;
        }
    }
}

Group File::createGroup(const std::string &path) {
    auto handle = H5Gcreate(self.handle, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    return Group(handle);
}

template<>
void
Group::write<double>(const std::string& dataSetName, const std::vector<dims_t>& dims, const double* data) {
    H5LTmake_dataset_double(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

}
}