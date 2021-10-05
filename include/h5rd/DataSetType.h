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
 * @file DataSetType.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include "Object.h"

namespace h5rd {
class DataSetType : public SubObject {
public:
    explicit DataSetType(handle_id hid, ParentFileRef parentFile);

    DataSetType(const DataSetType &rhs);

    DataSetType &operator=(const DataSetType &rhs);

    void close() override;

    ~DataSetType() override;
};

class VLENDataSetType : public DataSetType {
public:
    explicit VLENDataSetType(const DataSetType &other);
};

template<typename T>
class NativeDataSetType : public DataSetType {
public:
    explicit NativeDataSetType(ParentFileRef parentFile);

    using type = T;
};

template<typename T, unsigned int len>
class NativeArrayDataSetType : public DataSetType {
public:
    using type = typename std::remove_pointer<typename std::decay<T>::type>::type;

    explicit NativeArrayDataSetType(ParentFileRef parentFile);

    constexpr static unsigned int size = len;
private:
    NativeDataSetType<type> nativeType;
};

template<typename T, typename = typename std::enable_if<util::is_std_array<T>::value>::type>
class NativeStdArrayDataSetType : public DataSetType {
public:
    using type = typename T::value_type;

    NativeStdArrayDataSetType(ParentFileRef parentFile);

    constexpr static std::size_t size = std::tuple_size<T>::value;
private:
    NativeDataSetType<type> nativeType;
};

template<typename T>
class STDDataSetType : public DataSetType {
public:
    explicit STDDataSetType(ParentFileRef parentFile);

    using type = T;
};

template<typename T, unsigned int len>
class STDArrayDataSetType : public DataSetType {
public:
    using type = typename std::remove_pointer<typename std::decay<T>::type>::type;

    explicit STDArrayDataSetType(ParentFileRef parentFile);

    constexpr static unsigned int size = len;
private:
    STDDataSetType<type> stdType;
};

class NativeCompoundType : public DataSetType {
public:
    explicit NativeCompoundType(handle_id tid, ParentFileRef parentFile, std::vector<DataSetType> &&insertTypes);

private:
    std::vector<DataSetType> insertedTypes;
};

class STDCompoundType : public DataSetType {
public:
    explicit STDCompoundType(const NativeCompoundType &nativeType);
};

class NativeCompoundTypeBuilder {
public:
    explicit NativeCompoundTypeBuilder(std::size_t size, Object::ParentFileRef parentFile);

    NativeCompoundTypeBuilder &insert(const std::string &name, std::size_t offset, DataSetType &&type);

    template<typename T>
    NativeCompoundTypeBuilder &insert(const std::string &name, std::size_t offset);

    template<typename T, unsigned int size>
    NativeCompoundTypeBuilder &insertArray(const std::string &name, std::size_t offset);

    template<typename T, typename = typename std::enable_if<util::is_std_array<T>::value>::type>
    NativeCompoundTypeBuilder &insertStdArray(const std::string &name, std::size_t offset);

    NativeCompoundTypeBuilder &insertString(const std::string &name, std::size_t offset);

    NativeCompoundType build();

private:
    handle_id tid;
    Object::ParentFileRef _parentFile;
    std::vector<DataSetType> insertedTypes;
};

}

#include "detail/DataSetType_detail.h"
