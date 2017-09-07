/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * Header file for performance timers. Timers measure the time and write them to their target of type PerformanceData.
 * Multiple timers may write to the target at the same time, thus they obtain a lock on the target.
 *
 * Timers are created by PerformanceNodes. The nodes actually own the PerformanceData.
 * Nodes can create sub-nodes, a tree-like structure emerges.
 * A node may create many timers that write to its PerformanceData in a thread-safe manner, by locking the target.
 * On the other hand, creating sub-nodes is not designed for multi-threaded use.
 *
 * @file Timer.h
 * @brief Header file for performance timers.
 * @author clonker
 * @author chrisfroe
 * @date 13.07.16
 */

#pragma once

#include <memory>
#include <chrono>
#include <utility>

#include "logging.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

struct PerformanceData {
    using time = double;
    PerformanceData(time t, std::size_t c) : cumulativeTime(t), count(c) {}
    mutable time cumulativeTime = 0.;
    mutable std::size_t count = 0;
    mutable std::mutex mutex;
};

class Timer {
public:
    explicit Timer(const PerformanceData &target, bool measure);

    ~Timer() noexcept;

    Timer(const Timer &other) = delete;

    Timer(Timer &&other) = default;

    Timer &operator=(const Timer &other) = delete;

    Timer &operator=(Timer &&other) = delete;

private:
    bool measure;
    const PerformanceData &target;
    std::chrono::high_resolution_clock::time_point begin;
};

class PerformanceNode {
public:
    using performance_node_ref = std::unique_ptr<PerformanceNode>;
    using performance_mutex = std::mutex;
    using performance_lock = std::unique_lock<performance_mutex>;

    /**
     * Will hold a reference to itself as root.
     */
    PerformanceNode(const std::string &name, bool measure);

    /**
     * Will use root that was given.
     */
    PerformanceNode(const std::string &name, bool measure, const PerformanceNode &root);

    const PerformanceNode &subnode(const std::string &name) const;

    void clear();

    Timer timeit() const;

    const PerformanceData &data() const;

    const std::string &name() const;

    const PerformanceNode &direct_child(const std::string &name) const;

    /**
     * Get the child that corresponds to path with respect to *this. Path may contain '/' to get children of children ...
     * If path has a leading '/' the path will be resolved with respect to the root, not *this.
     * @param path
     * @return reference to desired node/child
     */
    const PerformanceNode &child(const std::string &path) const;

    const PerformanceNode &child(const std::vector<std::string> &labels) const;

    const std::size_t n_children() const;

    friend std::ostream& operator<<(std::ostream &os, const PerformanceNode &node) {
        os << node.describe();
        return os;
    }

    std::string describe(std::size_t level = 0) const;

    std::vector<std::string> keys() const;

private:
    std::vector<std::string> parse(const std::string &path) const;

    std::string validateName(const std::string &name) const;

    std::string _name;
    bool _measure;
    mutable std::vector<performance_node_ref> children;
    mutable performance_mutex childrenMutex;
    PerformanceData _data;
    std::reference_wrapper<const PerformanceNode> _root;
    static constexpr const char slash[] = "/";
};

NAMESPACE_END(util)
NAMESPACE_END(readdy)
