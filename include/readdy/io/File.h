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
 * @file File.h
 * @brief << brief description >>
 * @author clonker
 * @date 31.08.16
 */

#ifndef READDY_MAIN_FILE_H
#define READDY_MAIN_FILE_H

#include <string>
#include <vector>

namespace readdy {
namespace io {

template<typename T>
class DataSet;

class DataSetType;

class File;

class Group;

template<typename T>
class NativeDataSetType;

template<typename T>
class STDDataSetType;


class DataSetType {
    template<typename T>
    friend class DataSet;

public:
    using data_set_type_t = int;

    bool operator==(const DataSetType &) const;

    bool operator!=(const DataSetType &) const;

    using type = void;

    template<typename T>
    static data_set_type_t of_native(const std::vector<T> &) {
        return NativeDataSetType<T>().tid;
    }

    template<typename T>
    static data_set_type_t of_native(const T *const) {
        return NativeDataSetType<T>().tid;
    }

    template<typename T>
    static data_set_type_t of_std(const std::vector<T> &) {
        return STDDataSetType<T>().tid;
    }

    template<typename T>
    static data_set_type_t of_std(const T *const) {
        return STDDataSetType<T>().tid;
    }

    data_set_type_t getTID() const;

protected:
    DataSetType();

    data_set_type_t tid;
};

template<typename T>
class NativeDataSetType : public DataSetType {
public:
    NativeDataSetType();

    using type = T;
};

template<typename T>
class STDDataSetType : public DataSetType {
public:
    STDDataSetType();

    using type = T;
};

template<typename T>
class DataSet {
public:
    using handle_t = int;
    using dims_t = unsigned long long;
    static const dims_t UNLIMITED_DIMS;

    DataSet(const std::string &name, const Group &group, const std::vector<dims_t> &dims,
            const std::vector<dims_t> &maxDims);

    void append(const std::vector<dims_t> &dims, const T *data);

private:
    const std::vector<dims_t> &dims;
    const std::vector<dims_t> &maxDims;
    const Group &group;
    dims_t extensionDim;
    handle_t handle;
    handle_t memorySpace;
};

class Group {
    friend class File;

    template<typename T>
    friend class DataSet;

public:
    using handle_t = int;
    using dims_t = unsigned long long;

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<T> &data) {
        write(dataSetName, {data.size()}, data.data());
    }

    void write(const std::string &dataSetName, const std::string &string);

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<dims_t> &dims, const T *data);

    Group createGroup(const std::string &path);

protected:

    Group(const File &file) : Group(file, -1, "/") {};

    Group(const File &file, handle_t handle, const std::string &);

    handle_t handle;
    std::string path;
    const File &file;
};

class File {
    template<typename T>
    friend class DataSet;

public:

    enum class Action {
        CREATE, OPEN
    };

    enum class Flag {
        READ_ONLY = 0, READ_WRITE, OVERWRITE, FAIL_IF_EXISTS, CREATE_NON_EXISTING, DEFAULT /* = rw, create, truncate */
    };

    File(const std::string &path, const Action &action, const std::vector<Flag> &flag);

    File(const std::string &path, const Action &action, const Flag &flag = Flag::OVERWRITE);

    File(const File &) = delete;

    File &operator=(const File &) = delete;

    virtual ~File();

    void flush();

    void close();

    Group createGroup(const std::string &path);

    void write(const std::string &dataSetName, const std::string &data);

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<T> &data) {
        root.write(dataSetName, {data.size()}, data.data());
    }

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<Group::dims_t> &dims, const T *data) {
        root.write<T>(dataSetName, dims, data);
    }

private:
    std::string path_;
    Group root;
};
}
}
#endif //READDY_MAIN_FILE_H
