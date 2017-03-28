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
 * @file scoped_async.h
 * @brief << brief description >>
 * @author clonker
 * @date 28.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <future>
#include <readdy/common/macros.h>
#include <readdy/common/logging.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

class scoped_async {
    std::future<void> async_;
public:
    template<typename Function, typename... Args>
    explicit scoped_async(Function &&fun, Args &&... args)
            : async_(std::async(std::launch::async, std::forward<Function>(fun), std::forward<Args>(args)...)) {}

    ~scoped_async() {
        if(async_.valid()) async_.wait();
    }

    scoped_async(const scoped_async &) = delete;

    scoped_async &operator=(const scoped_async &) = delete;

    scoped_async(scoped_async &&) = default;

    scoped_async &operator=(scoped_async &&) = default;
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
