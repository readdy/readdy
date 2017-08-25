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
 * @file #.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 24.08.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/TopologyTypeRegistry.h>

unsigned short readdy::model::top::TopologyTypeRegistry::counter = 0;

void readdy::model::top::TopologyTypeRegistry::add_reaction(topology_type_type type,
                                                            const reactions::TopologyReaction &reaction) {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        it->second.internal_reactions.push_back(reaction);
    } else {
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }
}

void readdy::model::top::TopologyTypeRegistry::add_reaction(topology_type_type type,
                                                            reactions::TopologyReaction &&reaction) {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        it->second.internal_reactions.push_back(std::move(reaction));
    } else {
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }
}

const readdy::model::top::TopologyTypeInfo &
readdy::model::top::TopologyTypeRegistry::info_of(topology_type_type type) const {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        return it->second;
    }
    throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
}

const readdy::model::top::TopologyTypeInfo::internal_topology_reactions &
readdy::model::top::TopologyTypeRegistry::reactions_of(topology_type_type type) const {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        return it->second.internal_reactions;
    }
    log::warn("requested internal topology reactions of type {} which did not exist!", type);
    return defaultInfo.internal_reactions;
}

bool readdy::model::top::TopologyTypeRegistry::empty() {
    return _registry.empty();
}

bool readdy::model::top::TopologyTypeRegistry::contains_reactions() {
    return _contains_reactions;
}

void readdy::model::top::TopologyTypeRegistry::configure() {
    _contains_reactions = false;
    for(const auto& e : _registry) {
        _contains_reactions |= e.second.internal_reactions.empty();
    }
}

void readdy::model::top::TopologyTypeRegistry::debug_output() const {
    // todo
}

const std::string &readdy::model::top::TopologyTypeRegistry::name_of(readdy::topology_type_type type) const {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        return it->second.name;
    }
    throw std::invalid_argument(fmt::format("The requested type id {} did not exist.", type));
}

readdy::topology_type_type readdy::model::top::TopologyTypeRegistry::id_of(const std::string &name) const {
    using entry_type = topology_type_registry::value_type;
    auto it = std::find_if(_registry.begin(), _registry.end(), [&name](const entry_type& entry) {
        return entry.second.name == name;
    });
    if(it != _registry.end()) {
        return it->first;
    }
    throw std::invalid_argument(fmt::format("The requested type \"{}\" did not exist.", name));
}

readdy::topology_type_type readdy::model::top::TopologyTypeRegistry::add_type(const std::string &name,
                                                        const internal_topology_reactions &reactions) {
    using entry_type = topology_type_registry::value_type;
    auto it = std::find_if(_registry.begin(), _registry.end(), [&name](const entry_type& entry) {
        return entry.second.name == name;
    });
    if(it == _registry.end()) {
        const auto type = counter++;
        _registry[type] = {name, type, reactions};
        return type;
    }
    throw std::invalid_argument(fmt::format("The type {} already existed.", name));
}

void readdy::model::top::TopologyTypeRegistry::add_reaction(const std::string &type,
                                                            const reactions::TopologyReaction &reaction) {
    add_reaction(id_of(type), reaction);
}

void readdy::model::top::TopologyTypeRegistry::add_reaction(const std::string &type,
                                                            reactions::TopologyReaction &&reaction) {
    add_reaction(id_of(type), std::forward<reactions::TopologyReaction>(reaction));
}

const readdy::model::top::TopologyTypeRegistry::internal_topology_reactions &
readdy::model::top::TopologyTypeRegistry::reactions_of(const std::string &type) const {
    return reactions_of(id_of(type));
}
