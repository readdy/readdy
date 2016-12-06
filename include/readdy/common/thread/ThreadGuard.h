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
 * @file ThreadGuard.h
 * @brief << brief description >>
 * @author clonker
 * @date 01.08.16
 */

#ifndef READDY_CPUKERNEL_THREADGUARD_H
#define READDY_CPUKERNEL_THREADGUARD_H

#include <thread>

namespace readdy {
namespace util {
namespace thread {
class ThreadGuard {
    std::thread *t;
public:
    explicit ThreadGuard(std::thread *const t) : t(t) {}

    ~ThreadGuard() {
        if (t && t->joinable()) t->join();
    }

    ThreadGuard(const ThreadGuard &) = delete;

    ThreadGuard &operator=(const ThreadGuard &) = delete;

    ThreadGuard(ThreadGuard &&tg) : t(std::move(tg.t)) {
        tg.t = nullptr;
    };

    ThreadGuard &operator=(ThreadGuard &&tg) {
        t = std::move(tg.t);
        tg.t = nullptr;
        return *this;
    }
};
}
}
}
#endif //READDY_CPUKERNEL_THREADGUARD_H
