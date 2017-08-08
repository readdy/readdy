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
 * @file notification_barrier.h
 * @brief notofication_barrier header file
 * @author clonker
 * @date 24.11.16
 */
#pragma once

#include <mutex>
#include <condition_variable>
#include <atomic>
#include "../macros.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

/**
 * A notification barrier that lets all threads wait until another thread called ready().
 */
class notification_barrier {
public:
    /**
     * constructs a new notification barrier
     */
    explicit notification_barrier() = default;

    /**
     * waste time until ready() was called
     */
    void wait() const {
        std::unique_lock<std::mutex> lock(mutex);
        while(!done.load()) {
            cv.wait(lock, [this] { return done.load(); });
        }
    }

    /**
     * sets notification_barrier::done to true and notifies waiting threads
     */
    void ready() const {
        done.store(true);
        cv.notify_all();
    }

private:
    mutable std::mutex mutex;
    mutable std::condition_variable cv;
    mutable std::atomic<bool> done {false};
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
