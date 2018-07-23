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
