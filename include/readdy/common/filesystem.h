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
 * This header's purpose is mainly to replace boost's filesystem library. The dir_iterator gives access to files within
 * a directory using tinydir, the other methods are self-explanatory.
 *
 * @file filesystem.h
 * @brief Definition of the dir_iterator and some file system utilities.
 * @author clonker
 * @date 14.10.16
 */

#pragma once

#include <string>
#include <memory>
#include "macros.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(fs)

constexpr char separator =
#if READDY_WINDOWS
        '\\';
#else
        '/';
#endif

/**
 * Directory iterator class
 */
class dir_iterator {
public:
    /**
     * constructs a new directory iterator for a given path
     * @param path the path in which to look for directories
     */
    explicit dir_iterator(const std::string& path);
    /**
     * function to determine if there are more directories to come
     * @return true if there are more directories to be iterated through, false otherwise
     */
    bool has_next() const;
    /**
     * returns the name of the input directory without the leading path
     * @return the input directory name
     */
    std::string base_name() const;
    /**
     * advances the iterator by one
     * @return the current directory
     */
    std::string next();
    /**
     * destructor
     */
    virtual ~dir_iterator();

    /**
     * copying is not allowed
     */
    dir_iterator(const dir_iterator&) = delete;
    /**
     * copying is not allowed
     */
    dir_iterator& operator=(const dir_iterator&) = delete;
    /**
     * move constructor
     */
    dir_iterator(dir_iterator&&) noexcept;
    /**
     * move assign
     */
    dir_iterator& operator=(dir_iterator&&) noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

};

/**
 * returns the current system path
 * @return the current system path
 */
std::string current_path();

/**
 * checks if a path exists
 * @param path the path to check
 * @return true if it exists, otherwise false
 */
bool exists(const std::string& path);

/**
 * checks if the input path resolves to a file
 * @param path the path
 * @return true if it is indeed a file, otherwise false
 */
bool is_file(const std::string &path);

/**
 * checks if the input path resolves to a directory
 * @param path the path
 * @return true if it is indeed a directory, otherwise false
 */
bool is_directory(const std::string &path);

NAMESPACE_END(fs)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
