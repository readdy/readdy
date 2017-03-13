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
 * @file ScopedThread.h
 * @brief ScopedThread header file
 * @author clonker
 * @date 01.08.16
 */

#ifndef READDY_CPUKERNEL_SCOPEDTHREAD_H
#define READDY_CPUKERNEL_SCOPEDTHREAD_H

#include <thread>

namespace readdy {
namespace util {
namespace thread {

/**
 * scoped_thread implementation
 */
class scoped_thread {
    std::thread t;
public:
    /**
     * Creates a new scoped_thread based on a thread object
     * @param _t the reference thread
     */
    explicit scoped_thread(std::thread _t) : t(std::move(_t)) {
        if (!t.joinable()) throw std::logic_error("No thread!");
    }

    /**
     * joins the contained thread if its joinable
     */
    ~scoped_thread() {
        if (t.joinable()) {
            t.join();
        }
    }

    /**
     * moves another scoped_thread into this
     * @param rhs the other scoped_thread
     */
    scoped_thread(scoped_thread &&rhs) {
        t = std::move(rhs.t);
    }

    /**
     * moves another scoped_thread into this
     * @param rhs the other scoped_thread
     * @return myself
     */
    scoped_thread &operator=(scoped_thread &&rhs) {
        t = std::move(rhs.t);
        return *this;
    }

    /**
     * copying is not allowed
     */
    scoped_thread(const scoped_thread &) = delete;

    /**
     * copying is not allowed
     */
    scoped_thread &operator=(const scoped_thread &) = delete;
};
}
}
}
#endif //READDY_CPUKERNEL_SCOPEDTHREAD_H
