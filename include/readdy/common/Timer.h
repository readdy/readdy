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
#include "common.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

struct PerformanceData {
    /**
     * defines the type in which elapsed time is stored
     */
    using time = readdy::scalar;

    /**
     * Creates a new performance datum containing the cumulative time and the number of calls.
     * @param t the initial cumulative time
     * @param c the initial number of calls
     */
    PerformanceData(time t, std::size_t c) : _cumulativeTime(t), _count(c) {}

    /**
     * Record some elapsed time
     * @param elapsed the elapsed time
     */
    void record(time elapsed) const {
        std::unique_lock<std::mutex> lock(mutex);
        _cumulativeTime += elapsed;
        ++_count;
    }

    /**
     * the current value of the aggregated time in this datum
     * @return the cumulative time
     */
    time cumulativeTime() const {
        return _cumulativeTime;
    }

    /**
     * the number of calls that were made to record plus the initial number of counts given in the constructor
     * @return the number of calls
     */
    std::size_t count() const {
        return _count;
    }

    /**
     * clears this datum
     */
    void clear() const {
        std::unique_lock<std::mutex> lock(mutex);
        _cumulativeTime = 0.;
        _count = 0;
    }

private:
    mutable time _cumulativeTime = 0.;
    mutable std::size_t _count = 0;
    mutable std::mutex mutex;
};

class Timer {
public:
    /**
     * Creates a new timer object that will record the time elapsed between construction and destruction in a given
     * performance datum.
     *
     * @param target the performance datum to store the elapsed time in
     * @param measure if set to false, the measurements will not be stored
     */
    Timer(const PerformanceData &target, bool measure) : target(target), measure(measure) {
        if (measure) {
            begin = std::chrono::high_resolution_clock::now();
        }
    }

    /**
     * Destructor of the timer, recording the elapsed lifetime of this object into the given performance datum.
     */
    ~Timer() {
        if (measure) {
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            long elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
            auto elapsedSeconds =
                    static_cast<PerformanceData::time>(1e-6) * static_cast<PerformanceData::time>(elapsed);
            target.record(elapsedSeconds);
        }
    }

    /**
     * no copy
     */
    Timer(const Timer &) = delete;

    /**
     * no move
     */
    Timer(Timer &&) = default;

    /**
     * NO COPY ASSIGN
     */
    Timer &operator=(const Timer &) = delete;

    /**
     * NO MOVE ASSIGN
     */
    Timer &operator=(Timer &&) = delete;

private:
    bool measure;
    const PerformanceData &target;
    std::chrono::high_resolution_clock::time_point begin;
};

class PerformanceNode {
public:
    /**
     * a reference to a performance node
     */
    using performance_node_ref = std::unique_ptr<PerformanceNode>;
    /**
     * the mutex implementation to use
     */
    using performance_mutex = std::mutex;
    /**
     * the lock implementation to use
     */
    using performance_lock = std::unique_lock<performance_mutex>;

    /**
     * No-op performance node
     */
    PerformanceNode() : _name(""), _measure(false), _data(0, 0), _root(*this) {}

    /**
     * Will hold a reference to itself as root.
     */
    PerformanceNode(const std::string &name, bool measure)
            : _name(validateName(name)), _measure(measure), _data(0., 0), _root(*this) {}

    /**
     * Will use root that was given.
     */
    PerformanceNode(const std::string &name, bool measure, const PerformanceNode &root)
            : _name(validateName(name)), _measure(measure), _data(0., 0), _root(root) {}

    /**
     * Creates a subnode to this performance node
     * @param name subnode's name
     * @return the subnode
     */
    const PerformanceNode &subnode(const std::string &name) const;

    /**
     * clear all performance data in this node and the child nodes
     */
    void clear() {
        _data.clear();
        for (auto &c : children) {
            c->clear();
        }
    };

    /**
     * creates a RAII timer object belonging to this node
     * @return the timer object
     */
    Timer timeit() const {
        return Timer(_data, _measure);
    }

    /**
     * returns the performance datum of this
     * @return the performance datum
     */
    const PerformanceData &data() const {
        return _data;
    };

    /**
     * This node's name.
     * @return the node's name
     */
    const std::string &name() const {
        return _name;
    }

    /**
     * Yields the child node with the given name. Raises if no such node exists.
     * @param name the child node's name
     * @return the child node
     */
    const PerformanceNode &direct_child(const std::string &name) const {
        auto it = std::find_if(children.begin(), children.end(), [&name](const performance_node_ref &node) {
            return node->_name == name;
        });
        if (it == children.end()) {
            throw std::runtime_error(fmt::format("Child with name {} does not exist", name));
        }
        return **it;
    }

    /**
     * Get the child that corresponds to path with respect to *this. Path may contain '/' to get children of children ...
     * If path has a leading '/' the path will be resolved with respect to the root, not *this.
     * @param path
     * @return reference to desired node/child
     */
    const PerformanceNode &child(const std::string &path) const;

    /**
     * yields a child node
     * @param labels resolved path
     * @return the child node
     */
    const PerformanceNode &child(const std::vector<std::string> &labels) const;

    /**
     * returns the number of direct child nodes
     * @return the number of direct child nodes
     */
    const std::size_t n_children() const {
        return children.size();
    }

    /**
     * This is how we print this.
     * @param os the output stream
     * @param node a performance node instance
     * @return the output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const PerformanceNode &node) {
        os << node.describe();
        return os;
    }

    /**
     * Yields a string representation of this
     * @param level the indent level
     * @return a string representation
     */
    std::string describe(std::size_t level = 0) const;

    /**
     * Yields the child nodes' names.
     * @return the child nodes' names.
     */
    std::vector<std::string> keys() const {
        std::vector<std::string> k;
        std::transform(children.begin(), children.end(), std::back_inserter(k), [](const auto &child) -> std::string {
            return child->name();
        });
        return k;
    }

private:
    /**
     * Parse a path into a vector of strings, i.e., its components.
     * @param path the path
     * @return the path's components
     */
    std::vector<std::string> parse(const std::string &path) const;

    /**
     * validate a name
     * @param name the name
     * @return trimmed copy of the name if everything was o.k.
     */
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
