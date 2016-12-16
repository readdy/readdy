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
 * @file File.cpp
 * @brief Implementation file of File.h
 * @author clonker
 * @date 31.08.16
 */

#include "readdy/io/File.h"

#include <stdexcept>
#include <algorithm>

#include <hdf5_hl.h>
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


template<> const DataSetType::data_set_type_t STDDataSetType<short>::tid = H5Tcopy(H5T_STD_I16LE);
template<> const DataSetType::data_set_type_t STDDataSetType<unsigned short>::tid = H5Tcopy(H5T_STD_U16LE);
template<> const DataSetType::data_set_type_t STDDataSetType<int>::tid = H5Tcopy(H5T_STD_I32LE);
template<> const DataSetType::data_set_type_t STDDataSetType<unsigned int>::tid = H5Tcopy(H5T_STD_U32LE);
template<> const DataSetType::data_set_type_t STDDataSetType<long>::tid = H5Tcopy(H5T_STD_I64LE);
template<> const DataSetType::data_set_type_t STDDataSetType<unsigned long>::tid = H5Tcopy(H5T_STD_U64LE);
template<> const DataSetType::data_set_type_t STDDataSetType<long long>::tid = H5Tcopy(H5T_STD_I64LE);
template<> const DataSetType::data_set_type_t STDDataSetType<unsigned long long>::tid = H5Tcopy(H5T_STD_U64LE);
template<> const DataSetType::data_set_type_t STDDataSetType<float>::tid = H5Tcopy(H5T_IEEE_F32LE);
template<> const DataSetType::data_set_type_t STDDataSetType<double>::tid = H5Tcopy(H5T_IEEE_F64LE);
template<> const DataSetType::data_set_type_t STDDataSetType<std::string>::tid = H5Tcopy(H5T_C_S1);

template<> const DataSetType::data_set_type_t NativeDataSetType<short>::tid = H5Tcopy(H5T_NATIVE_SHORT);
template<> const DataSetType::data_set_type_t NativeDataSetType<unsigned short>::tid = H5Tcopy(H5T_NATIVE_USHORT);
template<> const DataSetType::data_set_type_t NativeDataSetType<int>::tid = H5Tcopy(H5T_NATIVE_INT);
template<> const DataSetType::data_set_type_t NativeDataSetType<unsigned int>::tid = H5Tcopy(H5T_NATIVE_UINT);
template<> const DataSetType::data_set_type_t NativeDataSetType<long>::tid = H5Tcopy(H5T_NATIVE_LONG);
template<> const DataSetType::data_set_type_t NativeDataSetType<unsigned long>::tid = H5Tcopy(H5T_NATIVE_ULONG);
template<> const DataSetType::data_set_type_t NativeDataSetType<long long>::tid = H5Tcopy(H5T_NATIVE_LLONG);
template<> const DataSetType::data_set_type_t NativeDataSetType<unsigned long long>::tid = H5Tcopy(H5T_NATIVE_ULLONG);
template<> const DataSetType::data_set_type_t NativeDataSetType<float>::tid = H5Tcopy(H5T_NATIVE_FLOAT);
template<> const DataSetType::data_set_type_t NativeDataSetType<double>::tid = H5Tcopy(H5T_NATIVE_DOUBLE);
template<> const DataSetType::data_set_type_t NativeDataSetType<bool>::tid = H5Tcopy(H5T_NATIVE_HBOOL);
template<> const DataSetType::data_set_type_t NativeDataSetType<std::string>::tid = H5Tcopy(H5T_C_S1);

Group::Group(const File &file, Group::handle_t handle, const std::string &path)
        : handle(handle), path(path), file(file) {}

void Group::write(const std::string &dataSetName, const std::string &string) {
    H5Tset_cset(STDDataSetType<std::string>::tid, H5T_CSET_UTF8);
    H5Tset_cset(NativeDataSetType<std::string>::tid, H5T_CSET_UTF8);
    H5LTmake_dataset_string(handle, dataSetName.c_str(), string.c_str());
}

Group Group::createGroup(const std::string &path) {
    auto handle = H5Gcreate(this->handle, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    return Group(file, handle, path);
}

Object::handle_t Group::getHandle() const {
    return handle;
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
        : path_(path), root{*this} {
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

const Group &File::getRootGroup() const {
    return root;
}

const unsigned long long Object::UNLIMITED_DIMS = H5S_UNLIMITED;

template<typename T>
void DataSet<T>::close() {
    if (memorySpace >= 0 && H5Sclose(memorySpace) < 0) {
        throw std::runtime_error("error on closing memory space!");
    } else {
        memorySpace = -1;
    }
    if (handle >= 0 && H5Dclose(handle) < 0) {
        throw std::runtime_error("error on closing data set!");
    } else {
        handle = -1;
    }
}

template<typename T>
DataSet<T>::~DataSet() {
    close();
}

template<typename T>
DataSet<T>::DataSet(const std::string &name, const Group &group, const std::vector<dims_t> &chunkSize,
                    const std::vector<dims_t> &maxDims) : maxDims(maxDims), group(group) {
    // validate and find extension dim
    {
        auto unlimited_it = std::find(maxDims.begin(), maxDims.end(), UNLIMITED_DIMS);
        bool containsUnlimited = unlimited_it != maxDims.end();
        if (!containsUnlimited) {
            throw std::runtime_error("needs to contain unlimited_dims in some dimension to be extensible");
        }
        extensionDim = static_cast<dims_t>(std::distance(maxDims.begin(), unlimited_it));
        log::console()->trace("found extension dim {}", extensionDim);
    }
    {
        // set up empty data set
        std::vector<dims_t> vec(maxDims.begin(), maxDims.end());
        vec[extensionDim] = 0;
        auto fileSpace = H5Screate_simple(static_cast<int>(maxDims.size()), vec.data(), maxDims.data());
        auto plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_chunk(plist, static_cast<int>(chunkSize.size()), chunkSize.data());
        handle = H5Dcreate(group.handle, name.c_str(), STDDataSetType<T>::tid, fileSpace, H5P_DEFAULT, plist, H5P_DEFAULT);
        memorySpace = -1;
        H5Pclose(plist);
        H5Pclose(fileSpace);
    }
}

template<typename T>
void DataSet<T>::append(const std::vector<dims_t> &dims, const T *data) {
    if(dims.size() != maxDims.size()) {
        log::console()->error("Tried to append data with ndims={} to set with ndims={}", dims.size(), maxDims.size());
        throw std::runtime_error("tried to append data with wrong dimensionality!");
    }
    if(memorySpace == -1) {
        memorySpace = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
    } else {
        H5Sset_extent_simple(memorySpace, static_cast<int>(dims.size()), dims.data(), nullptr);
    }
    std::vector<dims_t> currentExtent;
    std::vector<dims_t> offset;
    offset.resize(dims.size());
    {
        currentExtent.resize(dims.size());
        {
            auto fileSpace = H5Dget_space(handle);
            H5Sget_simple_extent_dims(fileSpace, currentExtent.data(), nullptr);
            H5Dclose(fileSpace);
        }
    }
    offset[extensionDim] = currentExtent[extensionDim];
    {
        std::vector<dims_t> newExtent(currentExtent);
        newExtent[extensionDim] += dims[extensionDim];
        H5Dset_extent(handle, newExtent.data());
    }
    auto fileSpace = H5Dget_space(handle);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.data(), nullptr, dims.data(), nullptr);
    H5Dwrite(handle, DataSetType::of_std(data), memorySpace, fileSpace, H5P_DEFAULT, data);
    H5Sclose(fileSpace);
}

DataSetType::DataSetType() {
}

}
}

template class readdy::io::DataSet<short>;
template class readdy::io::DataSet<int>;
template class readdy::io::DataSet<long>;
template class readdy::io::DataSet<float>;
template class readdy::io::DataSet<double>;