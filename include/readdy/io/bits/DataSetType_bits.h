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
 * @file DataSetType_bits.h
 * @brief << brief description >>
 * @author clonker
 * @date 04/01/2017
 * @copyright GNU Lesser General Public License v3.0
 */
#pragma once
#include <H5Tpublic.h>
#include <readdy/common/logging.h>
#include <readdy/common/traits.h>
#include "../DataSetType.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)

template<>
inline STDDataSetType<short>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_STD_I16LE)); }

template<>
inline STDDataSetType<unsigned short>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_STD_U16LE)); }

template<>
inline STDDataSetType<int>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_STD_I32LE)); }

template<>
inline STDDataSetType<unsigned int>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_STD_U32LE)); }

template<>
inline STDDataSetType<long>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_STD_I64LE)); }

template<>
inline STDDataSetType<unsigned long>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_STD_U64LE)); }

template<>
inline STDDataSetType<long long>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_STD_I64LE)); }

template<>
inline STDDataSetType<unsigned long long>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_STD_U64LE)); }

template<>
inline STDDataSetType<float>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_IEEE_F32LE)); }

template<>
inline STDDataSetType<double>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_IEEE_F64LE)); }

template<>
inline STDDataSetType<std::string>::STDDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_C_S1)); }

template<>
inline NativeDataSetType<short>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_SHORT)); }

template<>
inline NativeDataSetType<unsigned short>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_USHORT)); }

template<>
inline NativeDataSetType<int>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_INT)); }

template<>
inline NativeDataSetType<unsigned int>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_UINT)); }

template<>
inline NativeDataSetType<long>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_LONG)); }

template<>
inline NativeDataSetType<unsigned long>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_ULONG)); }

template<>
inline NativeDataSetType<long long>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_LLONG)); }

template<>
inline NativeDataSetType<unsigned long long>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_ULLONG)); }

template<>
inline NativeDataSetType<float>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_FLOAT)); }

template<>
inline NativeDataSetType<double>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_DOUBLE)); }

template<>
inline NativeDataSetType<bool>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_NATIVE_HBOOL)); }

template<>
inline NativeDataSetType<std::string>::NativeDataSetType() { tid = std::make_shared<DataTypeHandle>(H5Tcopy(H5T_C_S1)); }

template<typename T, typename enable>
inline NativeStdArrayDataSetType<T, enable>::NativeStdArrayDataSetType() {
    nativeType = NativeDataSetType<type>{};
    hsize_t dim[1] = {size};
    tid = std::make_shared<DataTypeHandle>(H5Tarray_create(nativeType.tid->tid, 1, dim));
}

template<typename T, unsigned int len>
inline NativeArrayDataSetType<T, len>::NativeArrayDataSetType() {
    nativeType = NativeDataSetType<type>{};
    hsize_t dim[1] = {len};
    tid = std::make_shared<DataTypeHandle>(H5Tarray_create(nativeType.tid->tid, 1, dim));
}

template<typename T, unsigned int len>
inline STDArrayDataSetType<T, len>::STDArrayDataSetType() {
    stdType = STDDataSetType<type>{};
    hsize_t dim[1] = {len};
    tid = std::make_shared<DataTypeHandle>(H5Tarray_create(stdType.tid->tid, 1, dim));
}

inline NativeCompoundType::NativeCompoundType(h5::data_set_type_t tid) {
    DataSetType::tid = std::make_shared<DataTypeHandle>(tid);
}

inline NativeCompoundTypeBuilder::NativeCompoundTypeBuilder(std::size_t size) {
    tid = H5Tcreate(H5T_COMPOUND, size);
}

inline NativeCompoundType NativeCompoundTypeBuilder::build() {
    return NativeCompoundType(tid);
}

inline NativeCompoundTypeBuilder& NativeCompoundTypeBuilder::insert(const std::string &name, std::size_t offset,
                                                 h5::data_set_type_t type) {
    if (H5Tinsert(tid, name.c_str(), offset, type) < 0) {
        H5Eprint(H5Eget_current_stack(), stderr);
        throw std::runtime_error(
                "error on inserting field " + name + " at offset " + std::to_string(offset) + " into compound type " +
                std::to_string(tid) + "!");
    }
    return *this;
}

template<typename T>
inline NativeCompoundTypeBuilder& NativeCompoundTypeBuilder::insert(const std::string& name, std::size_t offset) {
    NativeDataSetType<typename std::decay<T>::type> type;
    return insert(name, offset, type.tid->tid);
}

template<typename T, typename enable>
inline NativeCompoundTypeBuilder& NativeCompoundTypeBuilder::insertStdArray(const std::string &name,
                                                                            std::size_t offset) {
    io::NativeStdArrayDataSetType<T> type;
    return insert(name, offset, type.tid->tid);
}

template<typename T, unsigned int size>
inline NativeCompoundTypeBuilder& NativeCompoundTypeBuilder::insertArray(const std::string &name, std::size_t offset) {
    io::NativeArrayDataSetType<T, size> type;
    return insert(name, offset, type.tid->tid);
}

inline STDCompoundType::STDCompoundType(const NativeCompoundType &nativeType) {
    auto copy = H5Tcopy(nativeType.tid->tid);
    H5Tpack(copy);
    tid = std::make_shared<DataTypeHandle>(copy);
}

NAMESPACE_END(io)
NAMESPACE_END(readdy)
