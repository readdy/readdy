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

#include <readdy/model/topologies/TopologyRegistry.h>

namespace readdy {
namespace model {
namespace top {
unsigned short TopologyRegistry::counter = 0;

void TopologyRegistry::add_structural_reaction(topology_type_type type,
                                               const reactions::StructuralTopologyReaction &reaction) {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        it->second.structural_reactions.push_back(reaction);
    } else {
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }
}

void TopologyRegistry::add_structural_reaction(topology_type_type type,
                                               reactions::StructuralTopologyReaction &&reaction) {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        it->second.structural_reactions.push_back(std::move(reaction));
    } else {
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }
}

const TopologyTypeInfo &
TopologyRegistry::info_of(topology_type_type type) const {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        return it->second;
    }
    throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
}

const TopologyTypeInfo::structural_reaction_vector &
TopologyRegistry::reactions_of(topology_type_type type) const {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        return it->second.structural_reactions;
    }
    log::warn("requested internal topology reactions of type {} which did not exist!", type);
    return defaultInfo.structural_reactions;
}

bool TopologyRegistry::empty() {
    return _registry.empty();
}

bool TopologyRegistry::contains_structural_reactions() {
    return _contains_structural_reactions;
}

void TopologyRegistry::configure() {
    _topology_reaction_types.clear();
    _contains_structural_reactions = false;

    for(const auto& e : _registry) {
        _contains_structural_reactions |= e.second.structural_reactions.empty();
    }

    for(const auto& entry : _spatial_reactions) {
        _topology_reaction_types.emplace(std::get<0>(entry.first));
        _topology_reaction_types.emplace(std::get<1>(entry.first));
    }
}

void TopologyRegistry::debug_output() const {
    if(!_spatial_reactions.empty()) {
        log::debug(" - spatial topology reactions:");
        for(const auto& entry : _spatial_reactions) {
            for(const auto& reaction : entry.second) {
                auto d = fmt::format("SpatialTopologyReaction( ({} (id={}), {} (id={})) -> ({} (id={}), {} (id={})) , radius={}, rate={}",
                                     _type_registry.get().name_of(reaction.type1()), reaction.type1(),
                                     _type_registry.get().name_of(reaction.type2()), reaction.type2(),
                                     _type_registry.get().name_of(reaction.type_to1()), reaction.type_to1(),
                                     _type_registry.get().name_of(reaction.type_to2()), reaction.type_to2(),
                                     reaction.radius(), reaction.rate());
                log::debug("     * reaction {}", d);
            }
        }
        log::debug(" - structural topology reactions:");
        // todo

        log::debug(" - topology potential configuration:");
        // todo
    }
}

const std::string &TopologyRegistry::name_of(readdy::topology_type_type type) const {
    auto it = _registry.find(type);
    if(it != _registry.end()) {
        return it->second.name;
    }
    throw std::invalid_argument(fmt::format("The requested type id {} did not exist.", type));
}

readdy::topology_type_type TopologyRegistry::id_of(const std::string &name) const {
    using entry_type = type_registry::value_type;
    auto it = std::find_if(_registry.begin(), _registry.end(), [&name](const entry_type& entry) {
        return entry.second.name == name;
    });
    if(it != _registry.end()) {
        return it->first;
    }
    throw std::invalid_argument(fmt::format("The requested type \"{}\" did not exist.", name));
}

readdy::topology_type_type TopologyRegistry::add_type(const std::string &name,
                                                                              const structural_reactions &reactions) {
    using entry_type = type_registry::value_type;
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

void TopologyRegistry::add_structural_reaction(const std::string &type,
                                               const reactions::StructuralTopologyReaction &reaction) {
    add_structural_reaction(id_of(type), reaction);
}

void TopologyRegistry::add_structural_reaction(const std::string &type,
                                               reactions::StructuralTopologyReaction &&reaction) {
    add_structural_reaction(id_of(type), std::forward<reactions::StructuralTopologyReaction>(reaction));
}

const TopologyRegistry::structural_reactions &
TopologyRegistry::reactions_of(const std::string &type) const {
    return reactions_of(id_of(type));
}

TopologyRegistry::TopologyRegistry(
        const readdy::model::ParticleTypeRegistry &typeRegistry) : _type_registry(typeRegistry) {}



void TopologyRegistry::add_spatial_reaction(const std::string &name, const util::particle_type_pair &types,
                                            const util::particle_type_pair &types_to, scalar rate,
                                            scalar radius, reactions::STRMode mode) {
    spatial_reaction reaction(name, types, types_to, rate, radius, mode);
    validate_spatial_reaction(reaction);
    _spatial_reactions[types].push_back(std::move(reaction));
}

const TopologyRegistry::spatial_reaction_map &TopologyRegistry::spatial_reaction_registry() const {
    return _spatial_reactions;
}

void TopologyRegistry::add_spatial_reaction(const std::string &name, const std::string &typeFrom1,
                                            const std::string &typeFrom2, const std::string &typeTo1,
                                            const std::string &typeTo2, scalar rate, scalar radius,
                                            reactions::STRMode mode) {
    auto type_pair = std::make_tuple(_type_registry.get().id_of(typeFrom1), _type_registry.get().id_of(typeFrom2));
    auto type_pair_to = std::make_tuple(_type_registry.get().id_of(typeTo1), _type_registry.get().id_of(typeTo2));
    add_spatial_reaction(name, type_pair, type_pair_to, rate, radius, mode);
}

const TopologyRegistry::spatial_reactions &
TopologyRegistry::spatial_reactions_by_type(particle_type_type t1, particle_type_type t2) const {

    auto it = _spatial_reactions.find(std::tie(t1, t2));
    if(it != _spatial_reactions.end()) {
        return it->second;
    }
    return defaultTopologyReactions;
}

bool TopologyRegistry::is_spatial_reaction_type(particle_type_type type) const {
    return _topology_reaction_types.find(type) != _topology_reaction_types.end();
}

bool TopologyRegistry::is_spatial_reaction_type(const std::string &name) const {
    return is_spatial_reaction_type(_type_registry.get().id_of(name));
}

void TopologyRegistry::validate_spatial_reaction(const spatial_reaction &reaction) const {
    if(reaction.rate() <= 0) {
        throw std::invalid_argument("The rate of an external topology reaction (" + reaction.name() +
                                    ") should always be positive");
    }
    if(reaction.radius() <= 0) {
        throw std::invalid_argument("The radius of an external topology reaction("
                                    + reaction.name() + ") should always be positive");
    }
    auto info1 = _type_registry.get().info_of(reaction.type1());
    auto info2 = _type_registry.get().info_of(reaction.type2());
    auto infoTo1 = _type_registry.get().info_of(reaction.type_to1());
    auto infoTo2 = _type_registry.get().info_of(reaction.type_to2());

    if(info1.flavor != particleflavor::TOPOLOGY && info2.flavor != particleflavor::TOPOLOGY) {
        throw std::invalid_argument(
                fmt::format("At least one of the educt types ({} and {}) of topology reaction {} need "
                                    "to be topology flavored!", info1.name, info2.name, reaction.name())
        );
    }

    {
        using tr_entry = spatial_reaction_map::value_type;
        auto it = std::find_if(_spatial_reactions.begin(), _spatial_reactions.end(), [&reaction](const tr_entry &e) -> bool {
            return std::find_if(e.second.begin(), e.second.end(), [&reaction](const spatial_reaction& tr) {
                return tr.name() == reaction.name();
            }) != e.second.end();
        });
        if(it != _spatial_reactions.end()) {
            throw std::invalid_argument("An external topology reaction with the same name (" + reaction.name() +
                                        ") was already registered!");
        }
    }

    if(reaction.connect()) {
        if(infoTo1.flavor != particleflavor::TOPOLOGY || infoTo2.flavor != particleflavor::TOPOLOGY) {
            throw std::invalid_argument(
                    fmt::format("One or both of the target types ({} and {}) of topology reaction {} did not have "
                                        "the topology flavor and connect was set to true!",
                                infoTo1.name, infoTo2.name, reaction.name())
            );
        }
    } else {
        // flavors need to stay the same: topology particles must remain topology particles, same for normal particles
        if(info1.flavor != infoTo1.flavor) {
            throw std::invalid_argument(fmt::format("if connect is false, flavors need to remain same, which is not the case "
                                                            "for {} (flavor={}) -> {} (flavor={})", info1.name,
                                                    particleflavor::particle_flavor_to_str(info1.flavor), infoTo1.name,
                                                    particleflavor::particle_flavor_to_str(infoTo1.flavor)));
        }
        if(info2.flavor != infoTo2.flavor) {
            throw std::invalid_argument(fmt::format("if connect is false, flavors need to remain same, which is not the case "
                                                            "for {} (flavor={}) -> {} (flavor={})", info2.name,
                                                    particleflavor::particle_flavor_to_str(info2.flavor), infoTo2.name,
                                                    particleflavor::particle_flavor_to_str(infoTo2.flavor)));
        }
    }

}


api::PotentialConfiguration &TopologyRegistry::potential_configuration() {
    return potentialConfiguration_;
}

const api::PotentialConfiguration &TopologyRegistry::potential_configuration() const {
    return potentialConfiguration_;
}

void TopologyRegistry::configure_bond_potential(const std::string &type1, const std::string &type2,
                                                const api::Bond &bond) {
    potentialConfiguration_.pairPotentials[std::make_tuple(_type_registry.get().id_of(type1),
                                                           _type_registry.get().id_of(type2))].push_back(bond);
}

void TopologyRegistry::configure_angle_potential(const std::string &type1, const std::string &type2,
                                                 const std::string &type3, const api::Angle &angle) {
    potentialConfiguration_.anglePotentials[std::make_tuple(_type_registry.get().id_of(type1),
                                                            _type_registry.get().id_of(type2),
                                                            _type_registry.get().id_of(type3))].push_back(angle);
}

void TopologyRegistry::configure_torsion_potential(const std::string &type1, const std::string &type2,
                                                   const std::string &type3, const std::string &type4,
                                                   const api::TorsionAngle &torsionAngle) {
    potentialConfiguration_.torsionPotentials[std::make_tuple(_type_registry.get().id_of(type1),
                                                              _type_registry.get().id_of(type2),
                                                              _type_registry.get().id_of(type3),
                                                              _type_registry.get().id_of(type4))].push_back(torsionAngle);
}

}
}
}



