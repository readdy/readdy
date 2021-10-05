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
 * @file DataSetType_detail.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include <iostream>
#include <utility>
#include "../DataSetType.h"

inline h5rd::DataSetType::DataSetType(handle_id hid, ParentFileRef parentFile) : SubObject(std::move(parentFile)) {
    _hid = hid;
}

inline h5rd::DataSetType::DataSetType(const DataSetType &rhs) : SubObject(rhs._parentFile) {
    _hid = rhs._hid;
    _closed = rhs._closed;
    if (/*_parentFile && !_parentFile->closed() && */valid()) {
        if (H5Iinc_ref(_hid) < 0) {
            throw Exception("Error on increase HDF5 reference counter for data type");
        }
    }
}

inline h5rd::DataSetType &h5rd::DataSetType::operator=(const h5rd::DataSetType &rhs) {
    if (this != &rhs) {
        if (_hid != H5I_INVALID_HID) {
            if (H5Idec_ref(_hid) < 0) {
                throw Exception("Error in decrease of HDF5 reference counter for data type, copy assign");
            }
        }
        _hid = rhs._hid;
        _closed = rhs._closed;
        _parentFile = rhs._parentFile;
        if (/*_parentFile && !_parentFile->closed() && */valid()) {
            if (H5Iinc_ref(_hid) < 0) {
                throw Exception("Error in increase HDF5 reference counter for data type, copy assign");
            }
        }
    }
    return *this;
}

inline h5rd::DataSetType::~DataSetType() {
    try {
        close();
    } catch (const Exception &e) {
        std::cerr << "Unable to close hdf5 data type: " << e.what() << std::endl;
    }
}

inline void h5rd::DataSetType::close() {
    if (/*_parentFile && !_parentFile->closed() && */valid()) {
        if (H5Idec_ref(id()) < 0) {
            throw Exception("Error on decrease HDF5 reference counter for data type");
        }
    }
}

namespace h5rd {
inline VLENDataSetType::VLENDataSetType(const DataSetType &other) : DataSetType(H5Tvlen_create(other.id()),
                                                                                other.parentFile()) {}

template<>
inline STDDataSetType<unsigned char>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_I8LE),
                                                                                             std::move(parentFile)) {}

template<>
inline
STDDataSetType<short>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_I16LE),
                                                                              std::move(parentFile)) {}

template<>
inline STDDataSetType<unsigned short>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_U16LE),
                                                                                              std::move(parentFile)) {}

template<>
inline
STDDataSetType<int>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_I32LE),
                                                                            std::move(parentFile)) {}

template<>
inline STDDataSetType<unsigned int>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_U32LE),
                                                                                            std::move(parentFile)) {}

template<>
inline
STDDataSetType<long>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_I64LE),
                                                                             std::move(parentFile)) {}

template<>
inline STDDataSetType<unsigned long>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_U64LE),
                                                                                             std::move(parentFile)) {}

template<>
inline
STDDataSetType<long long>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_I64LE),
                                                                                  std::move(parentFile)) {}

template<>
inline
STDDataSetType<unsigned long long>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_STD_U64LE),
                                                                                           std::move(parentFile)) {}

template<>
inline
STDDataSetType<float>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_IEEE_F32LE),
                                                                              std::move(parentFile)) {}

template<>
inline
STDDataSetType<double>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_IEEE_F64LE),
                                                                               std::move(parentFile)) {}

template<>
inline
STDDataSetType<std::string>::STDDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_C_S1),
                                                                                    std::move(parentFile)) {}

template<>
inline
NativeDataSetType<unsigned char>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_UCHAR),
                                                                                            std::move(parentFile)) {}

template<>
inline NativeDataSetType<short>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_SHORT),
                                                                                           std::move(parentFile)) {}

template<>
inline
NativeDataSetType<unsigned short>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_USHORT),
                                                                                             std::move(parentFile)) {}

template<>
inline NativeDataSetType<int>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_INT),
                                                                                         std::move(parentFile)) {}

template<>
inline
NativeDataSetType<unsigned int>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_UINT),
                                                                                           std::move(parentFile)) {}

template<>
inline NativeDataSetType<long>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_LONG),
                                                                                          std::move(parentFile)) {}

template<>
inline
NativeDataSetType<unsigned long>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_ULONG),
                                                                                            std::move(parentFile)) {}

template<>
inline
NativeDataSetType<long long>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_LLONG),
                                                                                        std::move(parentFile)) {}

template<>
inline NativeDataSetType<unsigned long long>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(
        H5Tcopy(H5T_NATIVE_ULLONG), std::move(parentFile)) {}

template<>
inline NativeDataSetType<float>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_FLOAT),
                                                                                           std::move(parentFile)) {}

template<>
inline NativeDataSetType<double>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_DOUBLE),
                                                                                            std::move(parentFile)) {}

template<>
inline NativeDataSetType<bool>::NativeDataSetType(ParentFileRef parentFile) : DataSetType(H5Tcopy(H5T_NATIVE_HBOOL),
                                                                                          std::move(parentFile)) {}

template<>
inline NativeDataSetType<std::string>::NativeDataSetType(ParentFileRef parentFile)
    : DataSetType(H5Tcopy(H5T_C_S1), std::move(parentFile)) {}

template<typename T, typename enable>
inline NativeStdArrayDataSetType<T, enable>::NativeStdArrayDataSetType(ParentFileRef parentFile)
        : DataSetType(-1, parentFile), nativeType(NativeDataSetType < type > {parentFile}) {
    hsize_t dim[1] = {size};
    _hid = H5Tarray_create(nativeType.id(), 1, dim);
}

template<typename T, unsigned int len>
inline NativeArrayDataSetType<T, len>::NativeArrayDataSetType(ParentFileRef parentFile)
        : DataSetType(-1, parentFile), nativeType(NativeDataSetType < type > {parentFile}) {
    hsize_t dim[1] = {len};
    _hid = H5Tarray_create(nativeType.id(), 1, dim);
}

template<typename T, unsigned int len>
inline STDArrayDataSetType<T, len>::STDArrayDataSetType(ParentFileRef parentFile)
        : DataSetType(-1, parentFile), stdType(STDDataSetType<type>(parentFile)) {
    hsize_t dim[1] = {len};
    _hid = H5Tarray_create(stdType.id(), 1, dim);
}

inline NativeCompoundType::NativeCompoundType(handle_id tid, ParentFileRef parentFile,
                                              std::vector<DataSetType> &&insertTypes)
        : DataSetType(tid, std::move(parentFile)), insertedTypes(std::move(insertTypes)) {}

inline NativeCompoundTypeBuilder::NativeCompoundTypeBuilder(std::size_t size, Object::ParentFileRef parentFile) {
    tid = H5Tcreate(H5T_COMPOUND, size);
    _parentFile = std::move(parentFile);
}

inline NativeCompoundType NativeCompoundTypeBuilder::build() {
    return std::move(NativeCompoundType(tid, _parentFile, std::move(insertedTypes)));
}

inline NativeCompoundTypeBuilder &NativeCompoundTypeBuilder::insert(const std::string &name,
                                                                    std::size_t offset,
                                                                    DataSetType &&type) {
    insertedTypes.push_back(std::move(type));
    if (H5Tinsert(tid, name.c_str(), offset, insertedTypes.back().id()) < 0) {
        throw Exception(
                "error on inserting field " + name + " at offset " + std::to_string(offset) + " into compound type " +
                std::to_string(tid) + "!");
    }
    return *this;
}

template<typename T>
inline NativeCompoundTypeBuilder &
NativeCompoundTypeBuilder::insert(const std::string &name, std::size_t offset) {
    return insert(name, offset, NativeDataSetType<typename std::decay<T>::type>(_parentFile));
}

template<typename T, typename enable>
inline NativeCompoundTypeBuilder &NativeCompoundTypeBuilder::insertStdArray(const std::string &name,
                                                                            std::size_t offset) {
    return insert(name, offset, NativeStdArrayDataSetType<T>(_parentFile));
}

template<typename T, unsigned int size>
inline NativeCompoundTypeBuilder &
NativeCompoundTypeBuilder::insertArray(const std::string &name, std::size_t offset) {
    return insert(name, offset, NativeArrayDataSetType<T, size>(_parentFile));
}

inline NativeCompoundTypeBuilder &
NativeCompoundTypeBuilder::insertString(const std::string &name, std::size_t offset) {
    NativeDataSetType<std::string> t(_parentFile);
    H5Tset_cset(t.id(), H5T_CSET_UTF8);
    H5Tset_size(t.id(), H5T_VARIABLE);
    return insert(name, offset, std::move(t));
}

inline STDCompoundType::STDCompoundType(const NativeCompoundType &nativeType) : DataSetType(-1,
                                                                                            nativeType.parentFile()) {
    auto copy = H5Tcopy(nativeType.id());
    H5Tpack(copy);
    _hid = copy;
}
}
