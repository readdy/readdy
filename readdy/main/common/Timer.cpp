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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * @file Timer.cpp
 * @brief Implementation of timer classes
 * @author chrisfroe
 * @date 01.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/common/Timer.h>
#include <readdy/common/string.h>

namespace readdy {
namespace util {

constexpr const char PerformanceNode::slash[];

//todo slash is forbidden
PerformanceNode::PerformanceNode(const std::string &name, bool measure) : _name(name), _measure(measure), _data(0., 0) {}

//todo slash is forbidden
PerformanceNode &PerformanceNode::subnode(const std::string &name) {
    const auto &it = std::find_if(children.begin(), children.end(), [&name](const performance_node_ref &node) { return node->_name == name; });
    if (it == children.end()) {
        children.push_back(std::make_unique<PerformanceNode>(name, _measure));
        return *children.back();
    }
    return **it;
}

void PerformanceNode::clear() {
    _data.cumulativeTime = 0.;
    _data.count = 0;
    for (auto &c : children) {
        c->clear();
    }
}

Timer PerformanceNode::timeit() {
    return Timer(_data, _measure);
}

const PerformanceData &PerformanceNode::data() const {
    return _data;
}

const PerformanceNode &PerformanceNode::direct_child(const std::string &name) const {
    const auto &it = std::find_if(children.begin(), children.end(), [&name](const performance_node_ref &node) { return node->_name == name; });
    if (it == children.end()) {
        throw std::runtime_error(fmt::format("Child with name {} does not exist", name));
    }
    return **it;
}

const std::size_t PerformanceNode::n_children() const {
    return children.size();
}

std::string PerformanceNode::describe(const size_t level) const {
    std::string s;
    for (auto i = 0; i < level; ++i) {
        s += "\t";
    }
    s += fmt::format("{}: time {} s, count {}\n", _name, _data.cumulativeTime, _data.count);
    for (const auto &c : children) {
        s += c->describe(level + 1);
    }
    return s;
}

std::vector<std::string> PerformanceNode::parse(const std::string &path) const {
    std::vector<std::string> labels;
    std::string residualPath(path);
    auto slashPos = residualPath.find(slash);
    while (slashPos != residualPath.npos) {
        auto lhs = residualPath.substr(0, slashPos);
        auto rhs = residualPath.substr(slashPos + std::strlen(slash), residualPath.npos);
        util::str::trim(lhs);
        util::str::trim(rhs);
        residualPath = rhs;
        labels.push_back(lhs);
        slashPos = residualPath.find(slash);
    }
    labels.push_back(residualPath);
    return labels;
}

const PerformanceNode &PerformanceNode::child(const std::vector<std::string> &labels) const {
    if (labels.empty()) {
        throw std::invalid_argument("labels must not be empty");
    }
    std::vector<std::reference_wrapper<const PerformanceNode>> nodes;
    nodes.push_back(std::cref(direct_child(labels[0])));
    for (auto i = 1; i< labels.size(); ++i) {
        auto previousNode = nodes.back();
        auto nextNode = std::cref(previousNode.get().direct_child(labels[i]));
        nodes.push_back(nextNode);
    }
    return nodes.back();
}

const PerformanceNode &PerformanceNode::child(const std::string &path) const {
    if (path.find(slash)) {
        const auto &labels = parse(path);
        return child(labels);
    } else {
        return direct_child(path);
    }
}

const std::string &PerformanceNode::name() const {
    return _name;
}

Timer::Timer(PerformanceData &target, bool measure) : target(target), measure(measure) {
    if (measure) {
        begin = std::chrono::high_resolution_clock::now();
    }
}

Timer::~Timer() {
    if (measure) {
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        long elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        const auto elapsedSeconds = static_cast<PerformanceData::time>(1e-6) * static_cast<PerformanceData::time>(elapsed);
        std::unique_lock<std::mutex> lock(target.mutex);
        target.cumulativeTime += elapsedSeconds;
        target.count++;
    }
}

}
}