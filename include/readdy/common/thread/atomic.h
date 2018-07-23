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
 * This header file contains the definitions of:
 *     * copyable_atomic, a wrapped atomic variable that can be copied
 *
 * @file atomic.h
 * @brief definitions for everything that is related to std::atomic
 * @author clonker
 * @date 14.09.17
 * @copyright BSD-3
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