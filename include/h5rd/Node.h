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
 * @file Node.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include "common.h"

namespace h5rd {

template<typename Container>
class Node {
public:

    using FilterConfiguration = std::vector<Filter *>;

    Group createGroup(const std::string &path);

    std::vector<std::string> subgroups() const;

    std::vector<std::string> containedDataSets() const;

    Group getSubgroup(const std::string &name);

    std::shared_ptr<DataSet> getDataset(const std::string &name, const DataSetType &memoryType, const DataSetType &fileType);

    template<typename T>
    std::shared_ptr<DataSet> getDataset(const std::string &name);

    bool exists(const std::string &name) const;

    group_info info() const;

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<T> &data);

    void write(const std::string &dataSetName, const std::string &string);

    template<typename T>
    void write(const std::string &dataSetName, const dimensions &dims, const T *data);

    template<typename T>
    void read(const std::string &dataSetName, std::vector<T> &array, std::vector<hsize_t> stride = {});

    template<typename T>
    void readSelection(const std::string &dataSetName, std::vector<T> &array,
                       dimensions offsets = {}, dimensions stride = {}, dimensions count = {}, dimensions block = {});

    template<typename T>
    void read(const std::string &dataSetName, std::vector<T> &array, DataSetType *memoryType, DataSetType *fileType,
              std::vector<hsize_t> stride = {});

    template<typename T>
    void
    readSelection(const std::string &dataSetName, std::vector<T> &array, DataSetType *memoryType, DataSetType *fileType,
                  dimensions offsets = {}, dimensions stride = {}, dimensions count = {}, dimensions block = {});


    template<typename T>
    void readVLEN(const std::string &dataSetName, std::vector<std::vector<T>> &array);

    template<typename T>
    void readVLENSelection(const std::string &dataSetName, std::vector<std::vector<T>> &array,
                           dimensions offsets = {}, dimensions stride = {}, dimensions count = {},
                           dimensions block = {});

    template<typename T>
    void readVLEN(const std::string &dataSetName, std::vector<std::vector<T>> &array,
                  DataSetType *memoryType, DataSetType *fileType);

    template<typename T>
    void readVLENSelection(const std::string &dataSetName, std::vector<std::vector<T>> &array,
                           DataSetType *memoryType, DataSetType *fileType,
                           dimensions offsets = {}, dimensions stride = {}, dimensions count = {},
                           dimensions block = {});

    template<typename T>
    std::unique_ptr<DataSet> createDataSet(const std::string &name, const dimensions &chunkSize,
                                           const dimensions &maxDims, const FilterConfiguration &filters = {});

    std::unique_ptr<DataSet> createDataSet(const std::string &name, const dimensions &chunkSize,
                                           const dimensions &maxDims, const DataSetType &memoryType,
                                           const DataSetType &fileType, const FilterConfiguration &filters = {});

    template<typename T>
    std::unique_ptr<VLENDataSet> createVLENDataSet(const std::string &name, const dimensions &chunkSize,
                                                   const dimensions &maxDims, const FilterConfiguration &filters = {});

    std::unique_ptr<VLENDataSet> createVLENDataSet(const std::string &name, const dimensions &chunkSize,
                                                   const dimensions &maxDims, const DataSetType &memoryType,
                                                   const DataSetType &fileType,
                                                   const FilterConfiguration &filters = {});


private:

    std::vector<std::string> subElements(H5O_type_t type) const;

    Container *me();

    const Container *me() const;
};

}

#include "detail/Node_detail.h"
