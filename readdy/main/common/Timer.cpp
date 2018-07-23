/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * @file Timer.cpp
 * @brief Implementation of timer classes
 * @author chrisfroe
 * @author clonker
 * @date 01.09.17
 * @copyright GPL-3
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