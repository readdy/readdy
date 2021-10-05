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
 * @file Node_detail.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include <H5LTpublic.h>
#include <H5Tpublic.h>
#include <numeric>
#include <cmath>

#include "../Exception.h"
#include "../PropertyList.h"
#include "../Node.h"
#include "../Group.h"
#include "../DataSet.h"
#include "../common.h"

#include "String_utils.h"

template<typename Container>
inline bool h5rd::Node<Container>::exists(const std::string &name) const {
    return static_cast<bool>(H5Lexists(me()->id(), name.c_str(), H5P_DEFAULT));
}

template<typename Container>
inline h5rd::Group h5rd::Node<Container>::getSubgroup(const std::string &name) {
    auto id = me()->id();
    auto gid = H5Gopen(id, name.data(), H5P_DEFAULT);
    if (gid < 0) {
        throw Exception("Could not open subgroup " + name + "!");
    }
    Group group(name, me()->parentFile());
    group._hid = gid;
    return group;
}

template<typename Container>
inline h5rd::group_info h5rd::Node<Container>::info() const {
    auto id = me()->id();
    group_info info{};
    if (H5Gget_info(id, &info) < 0) {
        throw Exception("Failed to get group info!");
    }
    return info;
}

template<typename Container>
inline std::shared_ptr<h5rd::DataSet> h5rd::Node<Container>::getDataset(const std::string &name, const DataSetType &memoryType, const DataSetType &fileType) {
    auto dsId = H5Dopen(me()->id(), name.data(), H5P_DEFAULT);
    if (dsId < 0) {
        throw Exception("Could not open dataset " + name + "! (error code " + std::to_string(dsId) + ")");
    }
    auto ds = std::make_shared<DataSet>(me()->parentFile(), memoryType, fileType);
    ds->_hid = dsId;
    return ds;
}

template<typename Container>
template<typename T>
inline std::shared_ptr<h5rd::DataSet> h5rd::Node<Container>::getDataset(const std::string &name) {
    STDDataSetType<T> stdDST(me()->parentFile());
    NativeDataSetType<T> nDST(me()->parentFile());
    return getDataset(name, stdDST, nDST);
}

template<typename Container>
inline std::vector<std::string> h5rd::Node<Container>::subElements(H5O_type_t type) const {
    auto id = me()->id();
    std::vector<std::string> result;
    auto group_info = info();
    result.reserve(group_info.nlinks);
    for (std::size_t i = 0; i < group_info.nlinks; ++i) {
        H5O_info_t oinfo{};
        herr_t errorCode;
        # if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 12
        errorCode = H5Oget_info_by_idx(id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo, H5O_INFO_BASIC, H5P_DEFAULT);
        # else
        errorCode = H5Oget_info_by_idx(id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo, H5P_DEFAULT);
        # endif
        if (errorCode < 0) {
            throw std::runtime_error("failed to get info by idx, error code " + std::to_string(errorCode));
        }

        if (oinfo.type == type) {
            auto size = 1 + H5Lget_name_by_idx(id, ".", H5_INDEX_NAME, H5_ITER_INC, i, NULL, 0, H5P_DEFAULT);
            if (size < 0) {
                H5Eprint(H5Eget_current_stack(), stderr);
            }

            std::vector<char> c_string(static_cast<std::size_t>(size));
            H5Lget_name_by_idx(id, ".", H5_INDEX_NAME, H5_ITER_INC, i, c_string.data(), (std::size_t) size,
                               H5P_DEFAULT);
            std::string label(c_string.data());

            result.push_back(std::move(label));
        }
    }
    return result;
}

template<typename Container>
inline std::vector<std::string> h5rd::Node<Container>::containedDataSets() const {
    return subElements(H5O_TYPE_DATASET);
}

template<typename Container>
inline std::vector<std::string> h5rd::Node<Container>::subgroups() const {
    return subElements(H5O_TYPE_GROUP);
}

template<typename Container>
inline h5rd::Group h5rd::Node<Container>::createGroup(const std::string &path) {
    auto id = me()->id();
    Group group(path, me()->parentFile());
    if (util::groupExists(id, path)) {
        if ((group._hid = H5Gopen(id, path.c_str(), H5P_DEFAULT)) < 0) {
            throw Exception("Failed to open group (" + path + ")");
        }
        return group;
    } else {
        LinkCreatePropertyList plist(me()->parentFile());
        plist.set_create_intermediate_group();
        if ((group._hid = H5Gcreate(id, path.c_str(), plist.id(), H5P_DEFAULT, H5P_DEFAULT)) < 0) {
            throw Exception("Failed to create group with intermediates (" + path + ")");
        }
        return group;
    }
}

template<typename Container>
inline void h5rd::Node<Container>::write(const std::string &dataSetName, const std::string &string) {
    auto stdstr = STDDataSetType<std::string>(me()->parentFile());
    auto nativestr = NativeDataSetType<std::string>(me()->parentFile());
    H5Tset_cset(stdstr.id(), H5T_CSET_UTF8);
    H5Tset_cset(nativestr.id(), H5T_CSET_UTF8);

    if (H5LTmake_dataset_string(me()->id(), dataSetName.c_str(), string.c_str()) < 0) {
        throw Exception("there was a problem with writing \"" + string + "\" into a hdf5 file.");
    }
}

template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::write(const std::string &dataSetName, const std::vector<T> &data) {
    if (std::is_same<typename std::decay<T>::type, std::string>::value) {
        util::writeVector(me()->id(), dataSetName, data);
    } else {
        write(dataSetName, {data.size()}, data.data());
    }
}

namespace {
template<typename T>
inline void
help_write(h5rd::handle_id handle, const std::string &dsname, const h5rd::dimensions &dims, const T *data) {}

template<>
inline void
help_write(h5rd::handle_id handle, const std::string &dsname, const h5rd::dimensions &dims, const short *data) {
    H5LTmake_dataset_short(handle, dsname.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void
help_write(h5rd::handle_id handle, const std::string &dsname, const h5rd::dimensions &dims, const int *data) {
    H5LTmake_dataset_int(handle, dsname.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void
help_write(h5rd::handle_id handle, const std::string &dsname, const h5rd::dimensions &dims, const long *data) {
    H5LTmake_dataset_long(handle, dsname.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void
help_write(h5rd::handle_id handle, const std::string &dsname, const h5rd::dimensions &dims, const float *data) {
    H5LTmake_dataset_float(handle, dsname.data(), static_cast<int>(dims.size()), dims.data(), data);
}

template<>
inline void
help_write(h5rd::handle_id handle, const std::string &dsname, const h5rd::dimensions &dims, const double *data) {
    H5LTmake_dataset_double(handle, dsname.data(), static_cast<int>(dims.size()), dims.data(), data);
}
}


template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::write(const std::string &dataSetName, const dimensions &dims, const T *data) {
    using TT = typename std::decay<T>;
    help_write(me()->id(), dataSetName, dims, data);
}

template<typename Container>
template<typename T>
inline std::unique_ptr<h5rd::DataSet> h5rd::Node<Container>::createDataSet(
        const std::string &name, const dimensions &chunkSize, const dimensions &maxDims,
        const FilterConfiguration &filters) {
    return createDataSet(name, chunkSize, maxDims, STDDataSetType<T>(me()->parentFile()),
                         NativeDataSetType<T>(me()->parentFile()), filters);
}

template<typename Container>
inline std::unique_ptr<h5rd::DataSet>
h5rd::Node<Container>::createDataSet(const std::string &name, const h5rd::dimensions &chunkSize,
                                     const h5rd::dimensions &maxDims, const h5rd::DataSetType &memoryType,
                                     const h5rd::DataSetType &fileType, const FilterConfiguration &filters) {
    dimension extensionDim;
    {
        std::stringstream result;
        std::copy(maxDims.begin(), maxDims.end(), std::ostream_iterator<int>(result, ", "));
        // todo log::trace("creating data set with maxDims={}", result.str());
    }
    // validate and find extension dim
    {
        auto unlimited_it = std::find(maxDims.begin(), maxDims.end(), UNLIMITED_DIMS);
        bool containsUnlimited = unlimited_it != maxDims.end();
        if (!containsUnlimited) {
            throw std::runtime_error("needs to contain unlimited_dims in some dimension to be extensible");
        }
        extensionDim = static_cast<dimension>(std::distance(maxDims.begin(), unlimited_it));
        // todo log::trace("found extension dim {}", extensionDim);
    }
    handle_id handle;
    {
        // set up empty data set
        dimensions dims(maxDims.begin(), maxDims.end());
        dims[extensionDim] = 0;
        DataSpace fileSpace(me()->parentFile(), dims, maxDims);
        DataSetCreatePropertyList propertyList(me()->parentFile());
        propertyList.set_layout_chunked();
        propertyList.set_chunk(chunkSize);
        for (auto f : filters) {
            propertyList.activate_filter(f);
        }

        auto _hid = H5Dcreate(me()->id(), name.c_str(), fileType.id(),
                              fileSpace.id(), H5P_DEFAULT, propertyList.id(), H5P_DEFAULT);
        if (_hid < 0) {
            throw Exception("Error on creating data set " + std::to_string(_hid));
        } else {
            handle = _hid;
        }
    }
    auto ds = std::make_unique<DataSet>(me()->parentFile(), memoryType, fileType);
    ds->_hid = handle;
    ds->extensionDim() = extensionDim;
    return ds;
}

template<typename Container>
template<typename T>
inline std::unique_ptr<h5rd::VLENDataSet> h5rd::Node<Container>::createVLENDataSet(
        const std::string &name, const dimensions &chunkSize, const dimensions &maxDims,
        const FilterConfiguration &filters) {
    return createVLENDataSet(name, chunkSize, maxDims, STDDataSetType<T>(me()->parentFile()),
                             NativeDataSetType<T>(me()->parentFile()), filters);
}

template<typename Container>
inline std::unique_ptr<h5rd::VLENDataSet>
h5rd::Node<Container>::createVLENDataSet(const std::string &name, const h5rd::dimensions &chunkSize,
                                         const h5rd::dimensions &maxDims, const h5rd::DataSetType &memoryType,
                                         const h5rd::DataSetType &fileType, const FilterConfiguration &filters) {
    dimension extensionDim;
    VLENDataSetType vlenMemoryType(memoryType);
    VLENDataSetType vlenFileType(fileType);

    // validate and find extension dim
    {
        auto unlimited_it = std::find(maxDims.begin(), maxDims.end(), UNLIMITED_DIMS);
        bool containsUnlimited = unlimited_it != maxDims.end();
        if (!containsUnlimited) {
            throw std::runtime_error("needs to contain unlimited_dims in some dimension to be extensible");
        }
        extensionDim = static_cast<dimension>(std::distance(maxDims.begin(), unlimited_it));
    }
    handle_id hid;
    {
        // set up empty data set
        dimensions dims(maxDims.begin(), maxDims.end());
        dims[extensionDim] = 0;
        DataSpace fileSpace(me()->parentFile(), dims, maxDims);
        DataSetCreatePropertyList propertyList(me()->parentFile());
        propertyList.set_layout_chunked();
        propertyList.set_chunk(chunkSize);
        for (auto f : filters) {
            propertyList.activate_filter(f);
        }

        hid = H5Dcreate(me()->id(), name.c_str(), vlenFileType.id(),
                        fileSpace.id(), H5P_DEFAULT, propertyList.id(), H5P_DEFAULT);
        if (hid < 0) {
            throw Exception("Error on creating vlen data set " + std::to_string(hid));
        }
    }
    auto ds = std::make_unique<VLENDataSet>(me()->parentFile(), vlenMemoryType, vlenFileType);
    ds->extensionDim() = extensionDim;
    ds->_hid = hid;
    return ds;
}


template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::read(const std::string &dataSetName, std::vector<T> &array,
                                        std::vector<hsize_t> stride) {
    STDDataSetType<T> stdDST(me()->parentFile());
    NativeDataSetType<T> nDST(me()->parentFile());
    read(dataSetName, array, &stdDST, &nDST, stride);
}

template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::readSelection(const std::string &dataSetName, std::vector<T> &array,
                                                 dimensions offsets, dimensions stride, dimensions count,
                                                 dimensions block) {
    STDDataSetType<T> stdDST(me()->parentFile());
    NativeDataSetType<T> nDST(me()->parentFile());
    readSelection(dataSetName, array, &stdDST, &nDST, offsets, stride, count, block);
}

template<typename Container>
template<typename T>
inline void
h5rd::Node<Container>::readSelection(const std::string &dataSetName, std::vector<T> &array, DataSetType *memoryType,
                                     DataSetType *fileType,
                                     dimensions offsets, dimensions stride, dimensions count, dimensions block) {

    auto hid = H5Dopen(me()->id(), dataSetName.data(), H5P_DEFAULT);
    DataSpace fileSpace(me()->parentFile(), H5Dget_space(hid));

    if (count.empty()) {
        count.resize(fileSpace.ndim());
        const auto dims = fileSpace.dims();
        for (std::size_t d = 0; d < fileSpace.ndim(); ++d) {
            count[d] = static_cast<std::size_t>(std::ceil(dims[d] / static_cast<double>(stride[d])));
        }
    }

    DataSpace memorySpace(me()->parentFile(), count);

    const auto nElements = std::accumulate(count.begin(), count.end(), 1, std::multiplies<>());
    array.resize(static_cast<std::size_t>(nElements));

    H5Sselect_none(fileSpace.id());

    if (offsets.empty()) {
        offsets.resize(fileSpace.ndim());
    }

    auto status = H5Sselect_hyperslab(fileSpace.id(), H5S_SELECT_SET,
                                      offsets.data(),
                                      stride.empty() ? nullptr : stride.data(),
                                      count.data(),
                                      block.empty() ? nullptr : block.data());

    if (status < 0) {
        throw Exception("Failed selecting hyperslab for \"" + dataSetName + "\"!");
    }

    status = H5Dread(hid, memoryType->id(), memorySpace.id(), fileSpace.id(), H5P_DEFAULT, array.data());

    if (status < 0) {
        throw Exception("Failed reading \"" + dataSetName + "\"!");
    }

    H5Dclose(hid);
}

template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::read(const std::string &dataSetName, std::vector<T> &array, DataSetType *memoryType,
                                        DataSetType *fileType, std::vector<hsize_t> stride) {
    //blosc_compression::initialize();

    const auto n_array_dims = 1 + util::n_dims<T>::value;
    auto hid = H5Dopen(me()->id(), dataSetName.data(), H5P_DEFAULT);

    DataSpace fileSpace(me()->parentFile(), H5Dget_space(hid));

    const auto ndim = fileSpace.ndim();

    //if(ndim != n_array_dims) {
    //    log::error("wrong dimensionality: {} != {}", ndim, n_array_dims);
    //    throw std::invalid_argument("wrong dimensionality when attempting to read a data set");
    //}

    const auto dims = fileSpace.dims();

    if (!stride.empty()) {

        std::vector<h5rd::dimension> counts(dims.size());
        for (std::size_t d = 0; d < dims.size(); ++d) {
            counts[d] = static_cast<std::size_t>(std::ceil(dims[d] / static_cast<double>(stride[d])));
        }
        std::vector<hsize_t> offsets(dims.size(), 0);

        DataSpace memorySpace(me()->parentFile(), counts);

        const auto nElements = std::accumulate(counts.begin(), counts.end(), 1, std::multiplies<>());

        H5Sselect_none(fileSpace.id());
        array.resize(static_cast<std::size_t>(nElements));


        auto status = H5Sselect_hyperslab(fileSpace.id(), H5S_SELECT_SET, offsets.data(),
                                          stride.data(), counts.data(), nullptr);

        if (status < 0) {
            throw Exception("Failed selecting hyperslab for \"" + dataSetName + "\"!");
        }

        status = H5Dread(hid, memoryType->id(), memorySpace.id(), fileSpace.id(), H5P_DEFAULT, array.data());

        if (status < 0) {
            throw Exception("Failed reading \"" + dataSetName + "\"!");
        }

    } else {
        std::size_t required_length = 1;
        for (const auto dim : dims) {
            required_length *= dim;
        }
        H5Sselect_all(fileSpace.id());
        array.resize(required_length);

        auto result = H5Dread(hid, memoryType->id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data());

        if (result < 0) {
            throw Exception("Failed reading \"" + dataSetName + "\"!");
        }
    }

    H5Dclose(hid);
}

template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::readVLENSelection(const std::string &dataSetName, std::vector<std::vector<T>> &array,
                                                     h5rd::dimensions offsets, h5rd::dimensions stride,
                                                     h5rd::dimensions count, h5rd::dimensions block) {
    STDDataSetType<T> stdDST(me()->parentFile());
    NativeDataSetType<T> nDST(me()->parentFile());
    readVLENSelection(dataSetName, array, &stdDST, &nDST, offsets, stride, count, block);
}

template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::readVLENSelection(const std::string &dataSetName, std::vector<std::vector<T>> &array,
                                                     h5rd::DataSetType *memoryType, h5rd::DataSetType *fileType,
                                                     h5rd::dimensions offsets, h5rd::dimensions stride,
                                                     h5rd::dimensions count, h5rd::dimensions block) {
    auto hid = H5Dopen2(me()->id(), dataSetName.data(), H5P_DEFAULT);
    {
        DataSpace fileSpace(me()->parentFile(), H5Dget_space(hid));

        VLENDataSetType vlenMemoryType(*memoryType);
        VLENDataSetType vlenFileType(*fileType);

        if (count.empty()) {
            count.resize(fileSpace.ndim());
            const auto dims = fileSpace.dims();
            for (std::size_t d = 0; d < fileSpace.ndim(); ++d) {
                count[d] = static_cast<std::size_t>(std::ceil(dims[d] / static_cast<double>(stride[d])));
            }
        }

        DataSpace memorySpace(me()->parentFile(), count);

        const auto nElements = std::accumulate(count.begin(), count.end(), 1, std::multiplies<>());

        H5Sselect_none(fileSpace.id());

        if (offsets.empty()) {
            offsets.resize(fileSpace.ndim());
        }

        auto status = H5Sselect_hyperslab(fileSpace.id(), H5S_SELECT_SET,
                                          offsets.data(),
                                          stride.empty() ? nullptr : stride.data(),
                                          count.data(),
                                          block.empty() ? nullptr : block.data());
        if (status < 0) {
            throw Exception("Failed selecting VLEN hyperslab for \"" + dataSetName + "\"!");
        }

        {
            std::unique_ptr<hvl_t[]> rawData(new hvl_t[nElements]);

            status = H5Dread(hid, vlenMemoryType.id(), memorySpace.id(), fileSpace.id(), H5P_DEFAULT, rawData.get());
            //auto result = H5Dread(hid, vlenMemoryType.id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData.get());

            if (status < 0) {
                throw Exception("Failed reading \"" + dataSetName + "\"!");
            }

            array.reserve(static_cast<std::size_t>(nElements));
            for (auto it = rawData.get(); it != rawData.get() + nElements; ++it) {
                const auto &container = *it;
                auto ptr = static_cast<T *>(container.p);
                array.emplace_back(ptr, ptr + container.len);
            }
            H5Dvlen_reclaim(vlenMemoryType.id(), memorySpace.id(), H5P_DEFAULT, rawData.get());
        }
    }
    H5Dclose(hid);
}

template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::readVLEN(const std::string &dataSetName, std::vector<std::vector<T>> &array) {
    STDDataSetType<T> stdDST(me()->parentFile());
    NativeDataSetType<T> nDST(me()->parentFile());
    readVLEN(dataSetName, array, &stdDST, &nDST);
}

template<typename Container>
template<typename T>
inline void h5rd::Node<Container>::readVLEN(const std::string &dataSetName, std::vector<std::vector<T>> &array,
                                            h5rd::DataSetType *memoryType, h5rd::DataSetType *fileType) {

    VLENDataSetType vlenMemoryType(*memoryType);
    VLENDataSetType vlenFileType(*fileType);

    const auto n_array_dims = 1 + util::n_dims<T>::value;
    auto hid = H5Dopen2(me()->id(), dataSetName.data(), H5P_DEFAULT);

    DataSpace memorySpace(me()->parentFile(), H5Dget_space(hid));

    std::size_t required_length = 1;
    for (const auto dim : memorySpace.dims()) {
        required_length *= dim;
    }

    {
        std::unique_ptr<hvl_t[]> rawData(new hvl_t[required_length]);

        auto result = H5Dread(hid, vlenMemoryType.id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData.get());

        if (result < 0) {
            throw Exception("Failed reading \"" + dataSetName + "\"!");
        }

        array.resize(required_length);
        auto arrayIt = array.begin();
        for (auto it = rawData.get(); it != rawData.get() + required_length; ++it, ++arrayIt) {
            const auto &container = *it;
            auto ptr = static_cast<T *>(container.p);
            *arrayIt = std::vector<T>(ptr, ptr + container.len);
        }
        H5Dvlen_reclaim(vlenMemoryType.id(), memorySpace.id(), H5P_DEFAULT, rawData.get());
    }


    H5Dclose(hid);
}

template<typename Container>
inline Container *h5rd::Node<Container>::me() {
    return static_cast<Container *>(this);
}

template<typename Container>
inline const Container *h5rd::Node<Container>::me() const {
    return static_cast<const Container *>(this);
}
