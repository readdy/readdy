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
 * This header file contains the definition of joining_future, an RAII implementation of std::future.
 *
 * @file joining_future.h
 * @brief Definition of readdy::util::thread::joining_future
 * @author clonker
 * @date 13.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <future>
#include <readdy/common/macros.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

template<typename T>
class joining_future {
    std::future<T> _future;

public:
    /**
     * Generates a joining future out of a future
     * @param f the future
     */
    explicit joining_future(std::future<T> &&f) : _future(std::move(f)) {}

    /**
     * move is o.k.
     */
    joining_future(joining_future &&) noexcept = default;

    /**
     * move is o.k.
     */
    joining_future &operator=(joining_future &&) noexcept = default;

    /**
     * no copying
     */
    joining_future(const joining_future &) = delete;

    /**
     * no copying
     */
    joining_future &operator=(const joining_future &) = delete;

    /**
     * wait for the future if it is valid
     */
    ~joining_future() {
        if (_future.valid()) _future.wait();
    }
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
