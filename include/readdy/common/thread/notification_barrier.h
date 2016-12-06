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
 * @file notification_barrier.h
 * @brief << brief description >>
 * @author clonker
 * @date 24.11.16
 */

#ifndef READDY_MAIN_NOTIFICATION_BARRIER_H
#define READDY_MAIN_NOTIFICATION_BARRIER_H

#include <mutex>
#include <condition_variable>
#include <atomic>

namespace readdy {
namespace util {
namespace thread {

class notification_barrier {
public:
    explicit notification_barrier() { }

    /**
     * waste time until all threads reached the barrier
     */
    void wait() const {
        std::unique_lock<std::mutex> lock(mutex);
        while(!done.load()) {
            cv.wait(lock, [this] { return done.load(); });
        }
    }

    void ready() const {
        done.store(true);
        cv.notify_all();
    }

private:
    mutable std::mutex mutex;
    mutable std::condition_variable cv;
    mutable std::atomic<bool> done {false};
};

}
}
}
#endif //READDY_MAIN_NOTIFICATION_BARRIER_H
