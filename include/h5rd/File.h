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
 * @file File.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.17
 * @copyright BSD-3
 */

#pragma once

#include "Node.h"
#include "Object.h"

namespace h5rd {

class File : public Object, public Node<File>, public std::enable_shared_from_this<Object> {
public:

    enum class Action {
        CREATE, OPEN
    };

    enum class Flag {
        READ_ONLY = 0, READ_WRITE, OVERWRITE, FAIL_IF_EXISTS, CREATE_NON_EXISTING, DEFAULT /* = rw, create, truncate */
    };

    using Flags = std::vector<Flag>;

    static std::shared_ptr<File> open(const std::string &path, const Flag &flag);

    static std::shared_ptr<File> open(const std::string &path, const Flags &flags);

    static std::shared_ptr<File> create(const std::string &path, const Flag &flag);

    static std::shared_ptr<File> create(const std::string &path, const Flags &flags);

    File(const File &) = delete;

    File &operator=(const File &) = delete;

    File(File &&) = default;

    File &operator=(File &&) = default;

    ~File() override;

    void flush();

    void close() override;

    ParentFileRef ref() const;

protected:

    std::shared_ptr<Object> getptr();

    static void setUp(std::shared_ptr<File> file);

    File(const std::string &path, const Action &action, const Flags &flags);

    File(const std::string &path, const Action &action, const Flag &flag = Flag::OVERWRITE);

    std::string path;
    Action action;
    Flags flags;
};

}

#include "detail/File_detail.h"
