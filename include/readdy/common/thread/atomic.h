/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file atomic.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <atomic>
#include "../macros.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

template<typename T>
class copyable_atomic {
    std::atomic<T> _a;
public:
    copyable_atomic() : _a() {}

    explicit copyable_atomic(const std::atomic<T> &a) : _a(a.load()) {}

    copyable_atomic(const copyable_atomic &other) : _a(other._a.load()) {}

    copyable_atomic &operator=(const copyable_atomic &other) {
        _a.store(other._a.load());
        return *this;
    }

    const std::atomic<T> &operator*() const {
        return _a;
    }

    std::atomic<T> &operator*() {
        return _a;
    }

    const std::atomic<T> operator->() const {
        return &_a;
    }

    std::atomic<T> *operator->() {
        return &_a;
    }
};


NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)