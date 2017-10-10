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
 * This header file contains the definitions of:
 *     * copyable_atomic, a wrapped atomic variable that can be copied
 *
 * @file atomic.h
 * @brief definitions for everything that is related to std::atomic
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
    /**
     * the contained value
     */
    using value_type = T;

    /**
     * Creates a new copyable atomic of the specified type
     */
    copyable_atomic() : _a() {}

    /**
     * Instantiate this by an atomic. Not an atomic operation.
     * @param a the other atomic
     */
    explicit copyable_atomic(const std::atomic<T> &a) : _a(a.load()) {}

    /**
     * Copy constructor. Not an atomic CAS operation.
     * @param other the other copyable atomic
     */
    copyable_atomic(const copyable_atomic &other) : _a(other._a.load()) {}

    /**
     * Copy assign. Not an atomic CAS operation.
     * @param other the other copyable atomic
     * @return this
     */
    copyable_atomic &operator=(const copyable_atomic &other) {
        _a.store(other._a.load());
        return *this;
    }

    /**
     * Const dereference operator, yielding the underlying std atomic.
     * @return the underlying std atomic
     */
    const std::atomic<T> &operator*() const {
        return _a;
    }

    /**
     * Nonconst dereference operator, yielding the underlying std atomic.
     * @return the underlying std atomic
     */
    std::atomic<T> &operator*() {
        return _a;
    }

    /**
     * Const pointer-dereference operator, yielding a pointer to the underlying std atomic.
     * @return a pointer to the underlying std_atomic
     */
    const std::atomic<T> *operator->() const {
        return &_a;
    }

    /**
     * Nonconst pointer-dereference operator, yielding a pointer to the underlying std atomic.
     * @return a pointer to the underlying std_atomic
     */
    std::atomic<T> *operator->() {
        return &_a;
    }
};


NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)