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
 * @file semaphore.h
 * @brief header containing the semaphore class
 * @author clonker
 * @date 16.11.16
 */

#pragma once

#include <mutex>
#include <condition_variable>
#include "../macros.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

/**
 * counting semaphore implementation
 */
class semaphore {
public:
    /**
     * creates a new semaphore
     * @param count the count, defaults to 0
     */
    explicit semaphore(int count = 0) : count(count) {}

    /**
     * increase count
     */
    inline void signal() const {
        std::unique_lock<std::mutex> lock(mutex);
        count++;
        cv.notify_one();
    }

    /**
     * waste time while count is <= 0, then decrement
     */
    inline void wait() const {
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock, [this]() { return count > 0; });
        count--;
    }

    /**
     * waste time while count is <= 0
     */
    inline void wait_no_decrement() const {
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock, [this]() { return count > 0; });
    }

    semaphore(const semaphore&) = delete;
    semaphore(semaphore&&) = delete;
    semaphore& operator=(const semaphore&) = delete;
    semaphore& operator=(semaphore&&) = delete;

private:
    mutable std::mutex mutex;
    mutable std::condition_variable cv;
    mutable int count;
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
