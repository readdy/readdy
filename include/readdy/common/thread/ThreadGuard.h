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
