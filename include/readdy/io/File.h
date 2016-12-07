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

class Group {
public:
    using handle_t = int;
    using dims_t = unsigned long long;

    template<typename T>
    void write(const std::string& dataSetName, const std::vector<T> &data) {
        write(dataSetName, {data.size()}, data.data());
    }

    template<typename T>
    void write(const std::string& dataSetName, const std::vector<dims_t> &dims, const T* data);

protected:
    friend class File;
    Group() : Group(-1) {};
    Group(handle_t handle);
    handle_t handle;
};
class File {
public:

    enum class Action {
        CREATE, OPEN
    };

    enum class Flag {
        READ_ONLY, READ_WRITE, OVERWRITE, FAIL_IF_EXISTS, CREATE_NON_EXISTING
    };

    File(const std::string &path, const Action &action, const std::vector<Flag>& flag);
    File(const std::string &path, const Action &action, const Flag& flag);

    virtual ~File();

    void flush();

    void close();

    Group createGroup(const std::string& path);

    template<typename T>
    void write(const std::string& dataSetName, const std::vector<T> &data) {
        self.write(dataSetName, data.data(), {data.size()});
    }

    template<typename T>
    void write(const std::string& dataSetName, const std::vector<Group::dims_t> &dims, const T* data) {
        self.write<T>(dataSetName, dims, data);
    }

private:
    std::string path_;
    Group self;
};
}
}
#endif //READDY_MAIN_FILE_H
