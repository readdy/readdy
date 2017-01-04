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
#ifndef READDY_MAIN_DATASETTYPE_BITS_H
#define READDY_MAIN_DATASETTYPE_BITS_H

#include <H5Tpublic.h>
#include "../DataSetType.h"

namespace readdy {
namespace io {

template<>
inline STDDataSetType<short>::STDDataSetType() { tid = H5Tcopy(H5T_STD_I16LE); }

template<>
inline STDDataSetType<unsigned short>::STDDataSetType() { tid = H5Tcopy(H5T_STD_U16LE); }

template<>
inline STDDataSetType<int>::STDDataSetType() { tid = H5Tcopy(H5T_STD_I32LE); }

template<>
inline STDDataSetType<unsigned int>::STDDataSetType() { tid = H5Tcopy(H5T_STD_U32LE); }

template<>
inline STDDataSetType<long>::STDDataSetType() { tid = H5Tcopy(H5T_STD_I64LE); }

template<>
inline STDDataSetType<unsigned long>::STDDataSetType() { tid = H5Tcopy(H5T_STD_U64LE); }

template<>
inline STDDataSetType<long long>::STDDataSetType() { tid = H5Tcopy(H5T_STD_I64LE); }

template<>
inline STDDataSetType<unsigned long long>::STDDataSetType() { tid = H5Tcopy(H5T_STD_U64LE); }

template<>
inline STDDataSetType<float>::STDDataSetType() { tid = H5Tcopy(H5T_IEEE_F32LE); }

template<>
inline STDDataSetType<double>::STDDataSetType() { tid = H5Tcopy(H5T_IEEE_F64LE); }

template<>
inline STDDataSetType<std::string>::STDDataSetType() { tid = H5Tcopy(H5T_C_S1); }

template<>
inline NativeDataSetType<short>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_SHORT); }

template<>
inline NativeDataSetType<unsigned short>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_USHORT); }

template<>
inline NativeDataSetType<int>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_INT); }

template<>
inline NativeDataSetType<unsigned int>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_UINT); }

template<>
inline NativeDataSetType<long>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_LONG); }

template<>
inline NativeDataSetType<unsigned long>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_ULONG); }

template<>
inline NativeDataSetType<long long>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_LLONG); }

template<>
inline NativeDataSetType<unsigned long long>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_ULLONG); }

template<>
inline NativeDataSetType<float>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_FLOAT); }

template<>
inline NativeDataSetType<double>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_DOUBLE); }

template<>
inline NativeDataSetType<bool>::NativeDataSetType() { tid = H5Tcopy(H5T_NATIVE_HBOOL); }

template<>
inline NativeDataSetType<std::string>::NativeDataSetType() { tid = H5Tcopy(H5T_C_S1); }

}
}
#endif //READDY_MAIN_DATASETTYPE_BITS_H
