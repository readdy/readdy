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
 * @file barrier.h
 * @brief << brief description >>
 * @author clonker
 * @date 16.11.16
 */

#ifndef READDY_MAIN_BARRIER_H
#define READDY_MAIN_BARRIER_H

#include <mutex>
#include <condition_variable>

namespace readdy {
namespace util {
namespace thread {

class barrier {
public:
    explicit barrier(std::size_t count) : fallback(count), count(count), generation(0) { }

    /**
     * waste time until all threads reached the barrier
     */
    void wait() const {
        auto gen = generation;
        std::unique_lock<std::mutex> lock(mutex);
        if (--count == 0) {
            ++generation;
            count = fallback;
            cv.notify_all();
        } else {
            // wait while gen == generation, i.e., while count > 0
            cv.wait(lock, [this, gen] { return gen != generation; });
        }
    }

private:
    mutable std::mutex mutex;
    mutable std::condition_variable cv;
    mutable std::size_t fallback;
    mutable std::size_t count;
    mutable std::size_t generation;
};

}
}
}
#endif //READDY_MAIN_BARRIER_H
