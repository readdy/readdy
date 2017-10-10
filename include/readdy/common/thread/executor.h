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
 * This header file contains the definitions necessary for the executor, capable of executing tasks asynchronous with
 * an a priori selected strategy: Either pooled, with std::futures or std::async calls.
 *
 * @file executor.h
 * @brief definition of the readdy::util::thread::executor
 * @author clonker
 * @date 14.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <thread>
#include <vector>
#include <functional>

#include <readdy/common/macros.h>

#include "joining_future.h"
#include "scoped_thread.h"
#include "scoped_async.h"
#include "ctpl.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

NAMESPACE_BEGIN(executor_type)
struct pool {};
struct std_thread {};
struct std_async {};
NAMESPACE_END(executor_type)

class executor_base {
public:
    /**
     * a function that takes a thread-id as first argument and has no return value
     */
    using executable = std::function<void(std::size_t)>;

    /**
     * executes a vector of executables and then waits for all of them to finish
     * @param executables the executables
     */
    virtual void execute_and_wait(std::vector<executable> &&executables) const = 0;

    /**
     * Packs a function into a function that accepts the thread id and binds the other arguments.
     * @tparam F the function type
     * @tparam Args the argument types
     * @param f the function
     * @param args the arguments to be bound
     * @return an `executable`
     */
    template<typename F, typename... Args>
    executable pack(F &&f, Args &&... args) const {
        return std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Args>(args)...);
    }

};

template<typename executor_tag = readdy::util::thread::executor_type::pool>
class executor : public executor_base {
public:

    /**
     * Creates a new executor based on the executor_tag specified. If the tag is pool, a pool instance must be provided.
     * @param pool the pool if tag is pool, otherwise one can let it default to nullptr
     */
    explicit executor(ctpl::thread_pool *pool = nullptr) : pool(pool) {
        if (is_mode_pool() && pool == nullptr) {
            log::critical("selected execution type pool but did not pass a valid pool pointer to the executor!");
            throw std::invalid_argument(
                    "selected execution type pool but did not pass a valid pool pointer to the executor!");
        }
    };

    /**
     * default destructor
     */
    ~executor() = default;

    /**
     * move is o.k.
     */
    executor(executor &&) = default;

    /**
     * move is o.k.
     */
    executor &operator=(executor &&) = default;

    /**
     * no copying around
     */
    executor(const executor &) = delete;

    /**
     * no copying around
     */
    executor &operator=(const executor &) = delete;

    /**
     * checks if the executor mode is pool
     * @return true if the mode is pool, otherwise false
     */
    static constexpr bool is_mode_pool() {
        return std::is_same<executor_tag, executor_type::pool>::value;
    }

    /**
     * checks if the executor mode is std_thread
     * @return true if the mode is std_thread, otherwise false
     */
    static constexpr bool is_mode_thread() {
        return std::is_same<executor_tag, executor_type::std_thread>::value;
    }

    /**
     * checks if the executor mode is std_async
     * @return true if the mode is std_async, otherwise false
     */
    static constexpr bool is_mode_async() {
        return std::is_same<executor_tag, executor_type::std_async>::value;
    }

    /**
     * executes a bunch of executables, then waits for them to complete
     * @param executables the executables
     */
    void execute_and_wait(std::vector<executable> &&executables) const {
        execute_and_wait(std::forward<std::vector<executable>>(executables), executor_tag());
    };

private:
    mutable ctpl::thread_pool *pool {nullptr};
    mutable std::mutex mutex {};

    /**
     * std_async implementation of execute_and_wait
     * @param executables the executables
     */
    void execute_and_wait(std::vector<executable> &&executables, executor_type::std_async /*unused*/) const {
        std::vector<scoped_async> threads;
        threads.reserve(executables.size());

        std::size_t i = 0;
        for (auto &&executable : executables) {
            threads.emplace_back(std::move(executable), i++);
        }
    }

    /**
     * std_thread implementation of execute_and_wait
     * @param executables the executables
     */
    void execute_and_wait(std::vector<executable> &&executables, executor_type::std_thread /*unused*/) const {
        std::vector<scoped_thread> threads;
        threads.reserve(executables.size());

        std::size_t i = 0;
        for (auto &&executable : executables) {
            threads.emplace_back(std::move(executable), i++);
        }
    }

    /**
     * pool implementation of execute_and_wait
     * @param executables the executables
     */
    void execute_and_wait(std::vector<executable> &&executables, executor_type::pool /*unused*/) const {
        using future = joining_future<void>;

        std::unique_lock<std::mutex> lock(this->mutex);

        std::vector<future> futures;
        futures.reserve(executables.size());

        for (auto &&executable : executables) {
            futures.emplace_back(pool->push(std::move(executable)));
        }
    }
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
