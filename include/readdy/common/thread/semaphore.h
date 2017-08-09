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
