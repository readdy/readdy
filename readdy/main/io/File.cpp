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
 * @file File.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 31.08.16
 */

#include "readdy/io/File.h"

#include <stdexcept>

#include <hdf5_hl.h>
#include <algorithm>
#include <readdy/common/logging.h>

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
        case File::Flag::DEFAULT:
            return H5F_ACC_RDWR | H5F_ACC_CREAT | H5F_ACC_TRUNC;
    }
}

Group::Group(const File &file, Group::handle_t handle, const std::string &path)
        : handle(handle), path(path), file(file) {}

void Group::write(const std::string &dataSetName, const std::string &string) {
    H5LTmake_dataset_string(handle, dataSetName.c_str(), string.c_str());
}

Group Group::createGroup(const std::string &path) {
    auto handle = H5Gcreate(this->handle, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    return Group(file, handle, path);
}


Group File::createGroup(const std::string &path) {
    auto handle = H5Gcreate(root.handle, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    return Group(*this, handle, path);
}

template<>
void Group::write<short>(const std::string &dataSetName, const std::vector<dims_t> &dims, const short *data) {
    H5LTmake_dataset_short(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
void Group::write<int>(const std::string &dataSetName, const std::vector<dims_t> &dims, const int *data) {
    H5LTmake_dataset_int(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
void Group::write<long>(const std::string &dataSetName, const std::vector<dims_t> &dims, const long *data) {
    H5LTmake_dataset_long(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
void Group::write<float>(const std::string &dataSetName, const std::vector<dims_t> &dims, const float *data) {
    H5LTmake_dataset_float(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
void Group::write<double>(const std::string &dataSetName, const std::vector<dims_t> &dims, const double *data) {
    H5LTmake_dataset_double(handle, dataSetName.data(), static_cast<int>(dims.size()), dims.data(), data);
}

void File::close() {
    flush();
    if (root.handle >= 0 && H5Fclose(root.handle) < 0) {
        throw std::runtime_error("error when closing HDF5 file \"" + path_ + "\" with handle "
                                 + std::to_string(root.handle));
    } else {
        // we are now in closed state, set the handle to -1
        root.handle = -1;
    }
}

void File::flush() {
    if (root.handle >= 0 && H5Fflush(root.handle, H5F_SCOPE_LOCAL) < 0) {
        throw std::runtime_error("error when flushing HDF5 file \"" + path_ + "\" with handle "
                                 + std::to_string(root.handle));
    }
}

File::~File() {
    close();
}

File::File(const std::string &path, const Action &action, const Flag &flag) : File(path, action,
                                                                                   std::vector<Flag>{flag}) {}

File::File(const std::string &path, const File::Action &action, const std::vector<File::Flag> &flags)
        : path_(path), root({*this}) {
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

void File::write(const std::string &dataSetName, const std::string &data) {
    root.write(dataSetName, data);
}

template<typename T>
const typename DataSet<T>::dims_t DataSet<T>::UNLIMITED_DIMS = H5S_UNLIMITED;

template<typename T>
void setupDataSet(DataSet<T>& ds, const std::string &name, const hid_t group, const std::vector<typename DataSet<T>::dims_t> &dims,
                  const std::vector<typename DataSet<T>::dims_t> &maxDims, typename DataSet<T>::dims_t& extensionDim, typename DataSet<T>::handle_t& handle, typename DataSet<T>::handle_t& memorySpace) {
    // validate and find extension dim
    auto type = STDDataSetType<T>();
    {
        if (dims.size() != maxDims.size()) {
            throw std::runtime_error("dims and maxDims need to have same size!");
        }
        auto unlimited_it = std::find(maxDims.begin(), maxDims.end(), DataSet<T>::UNLIMITED_DIMS);
        bool containsUnlimited = unlimited_it != maxDims.end();
        if (!containsUnlimited) {
            throw std::runtime_error("needs to contain unlimited_dims in some dimension to be extensible");
        }
        extensionDim = static_cast<typename DataSet<T>::dims_t>(std::distance(maxDims.begin(), unlimited_it));
    }
    {
        // set up empty data set
        auto fileSpace = H5Screate_simple(static_cast<int>(dims.size()), std::vector<typename DataSet<T>::dims_t>(dims.size()).data(),
                                          maxDims.data());
        auto plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        // todo chunksize?
        handle = H5Dcreate(group, name.c_str(), type.getTID(), fileSpace, H5P_DEFAULT, plist, H5P_DEFAULT);
        memorySpace = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), NULL);
        H5Pclose(plist);
        H5Pclose(fileSpace);
    }
}

template<>
DataSet<short>::DataSet(const std::string &name, const Group &group, const std::vector<dims_t> &dims,
                         const std::vector<dims_t> &maxDims) : dims(dims), maxDims(maxDims), group(group) {
    setupDataSet(*this, name, group.handle, dims, maxDims, extensionDim, handle, memorySpace);
}

template<>
DataSet<int>::DataSet(const std::string &name, const Group &group, const std::vector<dims_t> &dims,
                         const std::vector<dims_t> &maxDims) : dims(dims), maxDims(maxDims), group(group) {
    setupDataSet(*this, name, group.handle, dims, maxDims, extensionDim, handle, memorySpace);
}

template<>
DataSet<long>::DataSet(const std::string &name, const Group &group, const std::vector<dims_t> &dims,
                         const std::vector<dims_t> &maxDims) : dims(dims), maxDims(maxDims), group(group) {
    setupDataSet(*this, name, group.handle, dims, maxDims, extensionDim, handle, memorySpace);
}

template<>
DataSet<float>::DataSet(const std::string &name, const Group &group, const std::vector<dims_t> &dims,
                         const std::vector<dims_t> &maxDims) : dims(dims), maxDims(maxDims), group(group) {
    setupDataSet(*this, name, group.handle, dims, maxDims, extensionDim, handle, memorySpace);
}

template<>
DataSet<double>::DataSet(const std::string &name, const Group &group, const std::vector<dims_t> &dims,
                 const std::vector<dims_t> &maxDims) : dims(dims), maxDims(maxDims), group(group) {
    setupDataSet(*this, name, group.handle, dims, maxDims, extensionDim, handle, memorySpace);
}

template<typename T>
hid_t
dataSetAppendPrepare(const std::vector<typename DataSet<T>::dims_t> &dims, typename DataSet<T>::handle_t memorySpace, typename DataSet<T>::handle_t handle,
                     typename DataSet<T>::dims_t extensionDim) {
    H5Sset_extent_simple(memorySpace, static_cast<int>(dims.size()), dims.data(), nullptr);
    std::vector<typename DataSet<T>::dims_t> currentExtent;
    {
        currentExtent.resize(dims.size());
        H5Sget_simple_extent_dims(handle, currentExtent.data(), nullptr);
    }
    {
        std::vector<typename DataSet<T>::dims_t> newExtent(currentExtent);
        newExtent[extensionDim] += dims[extensionDim];
        H5Dset_extent(handle, newExtent.data());
    }
    auto fileSpace = H5Dget_space(handle);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, currentExtent.data(), nullptr, dims.data(), nullptr);
    return fileSpace;
}

template<>
void DataSet<short>::append(const std::vector<dims_t> &dims, const short *data) {
    auto fileSpace = dataSetAppendPrepare<short>(dims, memorySpace, handle, extensionDim);
    H5Dwrite(handle, DataSetType::of_std(data), memorySpace, fileSpace, H5P_DEFAULT, data);
    H5Sclose(fileSpace);
}

template<>
void DataSet<int>::append(const std::vector<dims_t> &dims, const int *data) {
    auto fileSpace = dataSetAppendPrepare<int>(dims, memorySpace, handle, extensionDim);
    H5Dwrite(handle, DataSetType::of_std(data), memorySpace, fileSpace, H5P_DEFAULT, data);
    H5Sclose(fileSpace);
}

template<>
void DataSet<long>::append(const std::vector<dims_t> &dims, const long *data) {
    auto fileSpace = dataSetAppendPrepare<long>(dims, memorySpace, handle, extensionDim);
    H5Dwrite(handle, DataSetType::of_std(data), memorySpace, fileSpace, H5P_DEFAULT, data);
    H5Sclose(fileSpace);
}

template<>
void DataSet<float>::append(const std::vector<dims_t> &dims, const float *data) {
    auto fileSpace = dataSetAppendPrepare<float>(dims, memorySpace, handle, extensionDim);
    H5Dwrite(handle, DataSetType::of_std(data), memorySpace, fileSpace, H5P_DEFAULT, data);
    H5Sclose(fileSpace);
}

template<>
void DataSet<double>::append(const std::vector<dims_t> &dims, const double *data) {
    auto fileSpace = dataSetAppendPrepare<double>(dims, memorySpace, handle, extensionDim);
    H5Dwrite(handle, DataSetType::of_std(data), memorySpace, fileSpace, H5P_DEFAULT, data);
    H5Sclose(fileSpace);
}

DataSetType::DataSetType() : tid(-1) {
}

bool DataSetType::operator==(const DataSetType &rhs) const {
    return H5Tequal(tid, rhs.tid) > 0;
}

bool DataSetType::operator!=(const DataSetType &rhs) const {
    return !operator==(rhs);
}

DataSetType::data_set_type_t DataSetType::getTID() const {
    return tid;
}

template<>
STDDataSetType<short>::STDDataSetType() { tid = H5Tcopy(H5T_STD_I16LE); }

template<>
STDDataSetType<unsigned short>::STDDataSetType() { tid = H5Tcopy(H5T_STD_U16LE); }

template<>
STDDataSetType<int>::STDDataSetType() { tid = H5Tcopy(H5T_STD_I32LE); }

template<>
STDDataSetType<unsigned int>::STDDataSetType() { tid = H5Tcopy(H5T_STD_U32LE); }

template<>
STDDataSetType<long>::STDDataSetType() { tid = H5Tcopy(H5T_STD_I64LE); }

template<>
STDDataSetType<unsigned long>::STDDataSetType() { tid = H5Tcopy(H5T_STD_U64LE); }

template<>
STDDataSetType<long long>::STDDataSetType() { tid = H5Tcopy(H5T_STD_I64LE); }

template<>
STDDataSetType<unsigned long long>::STDDataSetType() { tid = H5Tcopy(H5T_STD_U64LE); }

template<>
STDDataSetType<float>::STDDataSetType() { tid = H5Tcopy(H5T_IEEE_F32LE); }

template<>
STDDataSetType<double>::STDDataSetType() { tid = H5Tcopy(H5T_IEEE_F64LE); }

template<>
STDDataSetType<std::string>::STDDataSetType() {
    tid = H5Tcopy(H5T_C_S1);
    H5Tset_cset(tid, H5T_CSET_UTF8);
}

template<>
NativeDataSetType<short>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_SHORT); }

template<>
NativeDataSetType<unsigned short>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_USHORT); }

template<>
NativeDataSetType<int>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_INT); }

template<>
NativeDataSetType<unsigned int>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_UINT); }

template<>
NativeDataSetType<long>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_LONG); }

template<>
NativeDataSetType<unsigned long>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_ULONG); }

template<>
NativeDataSetType<long long>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_LLONG); }

template<>
NativeDataSetType<unsigned long long>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_ULLONG); }

template<>
NativeDataSetType<float>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_FLOAT); }

template<>
NativeDataSetType<double>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_DOUBLE); }

template<>
NativeDataSetType<bool>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_HBOOL); }

template<>
NativeDataSetType<std::string>::NativeDataSetType() {
    tid = H5Tcopy(H5T_C_S1);
    H5Tset_cset(tid, H5T_CSET_UTF8);
}

}
}