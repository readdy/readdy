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
 * Header file containing the ThreadGuard class, which will - given a pointer to a thread - join that thread upon
 * destruction.
 *
 * @file ThreadGuard.h
 * @brief ThreadGuard header file
 * @author clonker
 * @date 01.08.16
 */

#pragma once

#include <thread>
#include "../macros.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

/**
 * Thread guard class that will, given a pointer to a thread, join that thread upon destruction.
 */
class ThreadGuard {
    std::thread *t;
public:
    /**
     * Constructs a new thread guard
     * @param t pointer to a thread object
     */
    explicit ThreadGuard(std::thread *const t) : t(t) {}

    /**
     * Joins the thread if it is not null and joinable
     */
    ~ThreadGuard() {
        if ((t != nullptr) && t->joinable()) t->join();
    }

    /**
     * Copying is not allowed
     */
    ThreadGuard(const ThreadGuard &) = delete;

    /**
     * Copying is not allowed
     */
    ThreadGuard &operator=(const ThreadGuard &) = delete;

    /**
     * Moves another ThreadGuard object into this one, setting the other's thread pointer to null in the process.
     * @param tg the other ThreadGuard object.
     */
    ThreadGuard(ThreadGuard &&tg) noexcept : t(tg.t) {
        tg.t = nullptr;
    };

    /**
     * See move constructor.
     * @param tg the other ThreadGuard object
     * @return myself
     */
    ThreadGuard &operator=(ThreadGuard &&tg) noexcept {
        t = tg.t;
        tg.t = nullptr;
        return *this;
    }
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
