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
 * @author clonker
 * @date 01.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/common/Timer.h>
#include <readdy/common/string.h>

namespace readdy {
namespace util {

constexpr const char PerformanceNode::slash[];

const PerformanceNode &PerformanceNode::subnode(const std::string &name) const {
    auto validatedName = validateName(name);
    performance_lock lock(childrenMutex);
    const auto &it = std::find_if(children.begin(), children.end(),
                                  [&validatedName](const performance_node_ref &node) {
                                      return node->_name == validatedName;
                                  });
    if (it == children.end()) {
        children.push_back(std::make_unique<PerformanceNode>(validatedName, _measure, _root.get()));
        return *children.back();
    }
    return **it;
}

std::string PerformanceNode::describe(const size_t level) const {
    std::string s;
    for (auto i = 0; i < level; ++i) {
        s += "\t";
    }
    s += fmt::format("{}: time {} s, count {}{}", _name, _data.cumulativeTime(), _data.count(),
                     readdy::util::str::newline);
    for (const auto &c : children) {
        s += c->describe(level + 1);
    }
    return s;
}

std::vector<std::string> PerformanceNode::parse(const std::string &path) const {
    std::vector<std::string> labels;
    std::string residualPath(path);
    auto slashPos = residualPath.find(slash);
    while (slashPos != std::string::npos) {
        auto lhs = residualPath.substr(0, slashPos);
        auto rhs = residualPath.substr(slashPos + std::strlen(slash), std::string::npos);
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
    std::reference_wrapper<const PerformanceNode> currentNode{direct_child(labels[0])};
    for (auto i = 1; i < labels.size(); ++i) {
        currentNode = std::cref(currentNode.get().direct_child(labels[i]));
    }
    return currentNode.get();
}

const PerformanceNode &PerformanceNode::child(const std::string &path) const {
    auto trimmedPath = util::str::trim_copy(path);
    if (trimmedPath.empty()) {
        return *this;
    }
    auto slashPos = trimmedPath.find(slash);
    if (slashPos != std::string::npos) {
        if (slashPos == 0) {
            auto rootPath = trimmedPath.substr(std::strlen(slash), std::string::npos);
            return _root.get().child(rootPath);
        } else {
            const auto &labels = parse(trimmedPath);
            return child(labels);
        }
    } else {
        return direct_child(trimmedPath);
    }
}

std::string PerformanceNode::validateName(const std::string &name) const {
    if (name.find(slash) != std::string::npos) {
        throw std::invalid_argument(fmt::format("name \"{}\" contains forbidden char {}", name, slash));
    }
    if (name.find(' ') == 0 || name.find(' ') == std::string::npos - 1) {
        throw std::invalid_argument(fmt::format("name \"{}\" contains leading/trailing whitespaces.", name));
    }
    return util::str::trim_copy(name);
}

}
}