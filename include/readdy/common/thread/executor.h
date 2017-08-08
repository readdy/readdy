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
 * @file executor.h
 * @brief << brief description >>
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
    using executable_t = std::function<void(std::size_t)>;

    virtual void execute_and_wait(std::vector<executable_t> &&executables) const = 0;

    template<typename F, typename... Args>
    executable_t pack(F &&f, Args &&... args) const {
        return std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Args>(args)...);
    }

};

template<typename executor_tag = readdy::util::thread::executor_type::pool>
class executor : public executor_base {
public:

    explicit executor(ctpl::thread_pool *pool = nullptr) : pool(pool) {
        if (is_mode_pool() && pool == nullptr) {
            log::critical("selected execution type pool but did not pass a valid pool pointer to the executor!");
            throw std::invalid_argument(
                    "selected execution type pool but did not pass a valid pool pointer to the executor!");
        }
    };

    ~executor() = default;

    executor(executor &&) = default;

    executor &operator=(executor &&) = default;

    executor(const executor &) = delete;

    executor &operator=(const executor &) = delete;

    static constexpr bool is_mode_pool() {
        return std::is_same<executor_tag, executor_type::pool>::value;
    }

    static constexpr bool is_mode_thread() {
        return std::is_same<executor_tag, executor_type::std_thread>::value;
    }

    static constexpr bool is_mode_async() {
        return std::is_same<executor_tag, executor_type::std_async>::value;
    }

    void execute_and_wait(std::vector<executable_t> &&executables) const {
        execute_and_wait(std::forward<std::vector<executable_t>>(executables), executor_tag());
    };

private:
    mutable ctpl::thread_pool *pool {nullptr};
    mutable std::mutex mutex {};

    void execute_and_wait(std::vector<executable_t> &&executables, executor_type::std_async /*unused*/) const {
        std::vector<scoped_async> threads;
        threads.reserve(executables.size());

        std::size_t i = 0;
        for (auto &&executable : executables) {
            threads.emplace_back(std::move(executable), i++);
        }
    }

    void execute_and_wait(std::vector<executable_t> &&executables, executor_type::std_thread /*unused*/) const {
        std::vector<scoped_thread> threads;
        threads.reserve(executables.size());

        std::size_t i = 0;
        for (auto &&executable : executables) {
            threads.emplace_back(std::move(executable), i++);
        }
    }

    void execute_and_wait(std::vector<executable_t> &&executables, executor_type::pool /*unused*/) const {
        using fut_t = joining_future<void>;

        std::unique_lock<std::mutex> lock(this->mutex);

        std::vector<fut_t> futures;
        futures.reserve(executables.size());

        for (auto &&executable : executables) {
            futures.emplace_back(pool->push(std::move(executable)));
        }
    }
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
