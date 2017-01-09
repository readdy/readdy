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
#include "Group.h"

namespace readdy {
namespace io {

class READDY_API File {
    template<typename T, bool VLEN>
    friend
    class DataSet;

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

    const Group &getRootGroup() const;

    void write(const std::string &dataSetName, const std::string &data);

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<T> &data) {
        root.write(dataSetName, {data.size()}, data.data());
    }

    template<typename T>
    void write(const std::string &dataSetName, const std::vector<h5::dims_t> &dims, const T *data) {
        root.write<T>(dataSetName, dims, data);
    }

private:
    std::string path_;
    Group root;
};

}
}

#include "bits/File_bits.h"

#endif //READDY_MAIN_FILE_H
