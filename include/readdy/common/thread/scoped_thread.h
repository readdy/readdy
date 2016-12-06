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
 * @file ScopedThread.h
 * @brief << brief description >>
 * @author clonker
 * @date 01.08.16
 */

#ifndef READDY_CPUKERNEL_SCOPEDTHREAD_H
#define READDY_CPUKERNEL_SCOPEDTHREAD_H

#include <thread>

namespace readdy {
namespace util {
namespace thread {

class scoped_thread {
    std::thread t;
public:
    explicit scoped_thread(std::thread _t) : t(std::move(_t)) {
        if (!t.joinable()) throw std::logic_error("No thread!");
    }

    ~scoped_thread() {
        if (t.joinable()) {
            t.join();
        }
    }

    scoped_thread(scoped_thread &&rhs) {
        t = std::move(rhs.t);
    }

    scoped_thread &operator=(scoped_thread &&rhs) {
        t = std::move(rhs.t);
        return *this;
    }

    scoped_thread(const scoped_thread &) = delete;

    scoped_thread &operator=(const scoped_thread &) = delete;
};
}
}
}
#endif //READDY_CPUKERNEL_SCOPEDTHREAD_H
