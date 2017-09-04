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
    time cumulativeTime = 0.;
    std::size_t count = 0;
    std::mutex mutex;
};

class Timer {
public:
    explicit Timer(PerformanceData &target, bool measure) : target(target), measure(measure) {
        if (measure) {
            begin = std::chrono::high_resolution_clock::now();
        }
    }

    ~Timer() noexcept {
        if (measure) {
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            long elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
            const auto elapsedSeconds = static_cast<PerformanceData::time>(1e-6) * static_cast<PerformanceData::time>(elapsed);
            std::unique_lock<std::mutex> lock(target.mutex);
            target.cumulativeTime += elapsedSeconds;
            target.count++;
        }
    }

    Timer(const Timer &other) = delete;

    Timer(Timer &&other) = default;

    Timer &operator=(const Timer &other) = delete;

    Timer &operator=(Timer &&other) = delete;

private:
    bool measure;
    PerformanceData &target;
    std::chrono::high_resolution_clock::time_point begin;
};

class PerformanceNode {
public:
    using performance_node_ref = std::unique_ptr<PerformanceNode>;
    PerformanceNode(const std::string &name, bool measure) : _name(name), _measure(measure), _data(0.,0) { }

    PerformanceNode &subnode(const std::string &name) {
        const auto& it = std::find_if(children.begin(), children.end(), [&name](const performance_node_ref& node){return node->_name == name;});
        if (it == children.end()) {
            children.push_back(std::make_unique<PerformanceNode>(name, _measure));
            return *children.back();
        }
        return **it;
    }

    void clear() {
        _data.cumulativeTime = 0.;
        _data.count = 0;
        for (auto &c : children) {
            c->clear();
        }
    }

    Timer timeit() {
        return Timer(_data, _measure);
    }

    const PerformanceData &data() const {
        return _data;
    }

    const PerformanceNode &child(const std::string &name) const {
        const auto& it = std::find_if(children.begin(), children.end(), [&name](const performance_node_ref& node){return node->_name == name;});
        if (it == children.end()) {
            throw std::runtime_error(fmt::format("Child with name {} does not exist", name));
        }
        return **it;
    }

    const std::size_t n_children() const {
        return children.size();
    }

private:
    std::string _name;
    bool _measure;
    std::vector<performance_node_ref> children;
    PerformanceData _data;
};

NAMESPACE_END(util)
NAMESPACE_END(readdy)
