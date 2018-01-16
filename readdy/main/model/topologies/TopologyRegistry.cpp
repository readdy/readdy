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
#include <readdy/model/Utils.h>

namespace readdy {
namespace model {
namespace top {
topology_type_type TopologyRegistry::counter = 0;

void TopologyRegistry::addStructuralReaction(topology_type_type type,
                                             const reactions::StructuralTopologyReaction &reaction) {
    auto it = _registry.find(type);
    if (it != _registry.end()) {
        it->second.structural_reactions.push_back(reaction);
    } else {
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }
}

void TopologyRegistry::addStructuralReaction(topology_type_type type,
                                             reactions::StructuralTopologyReaction &&reaction) {
    auto it = _registry.find(type);
    if (it != _registry.end()) {
        it->second.structural_reactions.push_back(std::move(reaction));
    } else {
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }
}

void TopologyRegistry::configure() {
    _topologyReactionTypes.clear();
    _containsStructuralReactions = false;

    for (const auto &e : _registry) {
        _containsStructuralReactions |= e.second.structural_reactions.empty();
    }

    for (const auto &entry : _spatialReactions) {
        _topologyReactionTypes.emplace(std::get<0>(entry.first));
        _topologyReactionTypes.emplace(std::get<2>(entry.first));
    }
}

std::string TopologyRegistry::describe() const {
    namespace rus = readdy::util::str;
    std::string description;
    if (!_spatialReactions.empty()) {
        description += fmt::format(" - spatial topology reactions:{}", rus::newline);
        for (const auto &entry : _spatialReactions) {
            for (const auto &reaction : entry.second) {
                auto d = fmt::format(
                        "SpatialTopologyReaction( ({} (id={}), {} (id={})) -> ({} (id={}), {} (id={})) , radius={}, rate={}",
                        _typeRegistry.get().nameOf(reaction.type1()), reaction.type1(),
                        _typeRegistry.get().nameOf(reaction.type2()), reaction.type2(),
                        _typeRegistry.get().nameOf(reaction.type_to1()), reaction.type_to1(),
                        _typeRegistry.get().nameOf(reaction.type_to2()), reaction.type_to2(),
                        reaction.radius(), reaction.rate());
                description += fmt::format("     * reaction {}{}", d, rus::newline);
            }
        }
    }
    description += fmt::format(" - topology types:{}", rus::newline);
    for (const auto &entry : _registry) {
        description += fmt::format("     * topology type \"{}\" with id {} and {} structural reactions{}",
                                   entry.second.name, entry.second.type, entry.second.structural_reactions.size(),
                                   rus::newline);
    }

    description += fmt::format(" - structural topology reactions:{}", rus::newline);
    for (const auto &entry : _registry) {
        description += fmt::format("     - topology type \"{}\" with {} structural reactions:{}",
                                   entry.second.name, entry.second.structural_reactions.size(), rus::newline);
        for (const auto &r : entry.second.structural_reactions) {
            description += fmt::format("         * reaction with roll_back = {} and create child tops = {}{}",
                                       r.rolls_back_if_invalid(), r.creates_child_topologies_after_reaction(), rus::newline);
        }
    }

    description += fmt::format(" - topology potential configuration:{}", rus::newline);
    description += fmt::format("     - bonds ({}):{}", _potentialConfiguration.pairPotentials.size(), rus::newline);
    for (const auto &entry : _potentialConfiguration.pairPotentials) {
        description += fmt::format("         - Bonds for particle types {} and {}:{}",
                                   _typeRegistry.get().nameOf(std::get<0>(entry.first)),
                                   _typeRegistry.get().nameOf(std::get<1>(entry.first)), rus::newline);
        auto bondToStr = [](const api::Bond &bond) -> std::string {
            switch (bond.type) {
                case api::BondType::HARMONIC:
                    return "Harmonic";
            }
        };
        for (const auto &bond : entry.second) {
            description += fmt::format("             * {} bond with force constant {} and length {}{}", bondToStr(bond),
                                       bond.forceConstant, bond.length, rus::newline);
        }
    }
    description += fmt::format("     - angles ({}):{}", _potentialConfiguration.anglePotentials.size(), rus::newline);
    for (const auto &entry : _potentialConfiguration.anglePotentials) {
        auto angleToStr = [](const api::Angle &angle) -> std::string {
            switch (angle.type) {
                case api::AngleType::HARMONIC:
                    return "Harmonic";
            }
        };
        for (const auto &angle : entry.second) {
            description += fmt::format("             * {} angle with force constant {} and equilibrium angle {}{}",
                                       angleToStr(angle), angle.forceConstant, angle.equilibriumAngle, rus::newline);
        }
    }
    description += fmt::format("     - torsions ({}):{}", _potentialConfiguration.torsionPotentials.size(), rus::newline);
    for (const auto &entry : _potentialConfiguration.torsionPotentials) {
        auto torsionToStr = [](const api::TorsionAngle &torsion) -> std::string {
            switch (torsion.type) {
                case api::TorsionType::COS_DIHEDRAL:
                    return "Cosine-Dihedral";
            }
        };
        for (const auto &dih : entry.second) {
            description += fmt::format(
                    "             * {} with force constant {}, equilibrium angle {} and multiplicity {}{}",
                    torsionToStr(dih), dih.forceConstant, dih.phi_0, dih.multiplicity, rus::newline);
        }
    }
    return description;
}

readdy::topology_type_type TopologyRegistry::addType(const std::string &name, const structural_reactions &reactions) {
    util::validateTypeName(name);
    using entry_type = type_registry::value_type;
    auto it = std::find_if(_registry.begin(), _registry.end(), [&name](const entry_type &entry) {
        return entry.second.name == name;
    });
    if (it == _registry.end()) {
        const auto type = counter++;
        _registry[type] = {name, type, reactions};
        return type;
    }
    throw std::invalid_argument(fmt::format("The type {} already existed.", name));
}

TopologyRegistry::TopologyRegistry(
        const readdy::model::ParticleTypeRegistry &typeRegistry) : _typeRegistry(typeRegistry) {}


void TopologyRegistry::addSpatialReaction(const std::string &name, const readdy::util::particle_type_pair &types,
                                          const topology_type_pair &topology_types,
                                          const readdy::util::particle_type_pair &types_to,
                                          const topology_type_pair &topology_types_to,
                                          scalar rate, scalar radius, reactions::STRMode mode) {
    spatial_reaction reaction(name, types, topology_types, types_to, topology_types_to, rate, radius, mode);
    addSpatialReaction(std::move(reaction));
}

void TopologyRegistry::addSpatialReaction(const std::string &descriptor, scalar rate, scalar radius) {
    reactions::STRParser parser(*this);
    addSpatialReaction(parser.parse(descriptor, rate, radius));
}

void TopologyRegistry::addSpatialReaction(reactions::SpatialTopologyReaction &&reaction) {
    validateSpatialReaction(reaction);
    auto key = std::make_tuple(reaction.type1(), reaction.top_type1(), reaction.type2(), reaction.top_type2());
    _spatialReactions[key].push_back(std::move(reaction));
}

void TopologyRegistry::addSpatialReaction(const std::string &name, const std::string &typeFrom1,
                                          const std::string &typeFrom2, const std::string &topologyTypeFrom1,
                                          const std::string &topologyTypeFrom2, const std::string &typeTo1,
                                          const std::string &typeTo2, const std::string &topologyTypeTo1,
                                          const std::string &topologyTypeTo2, scalar rate, scalar radius,
                                          reactions::STRMode mode) {
    auto type_pair = std::make_tuple(_typeRegistry.get().idOf(typeFrom1), _typeRegistry.get().idOf(typeFrom2));
    auto type_pair_to = std::make_tuple(_typeRegistry.get().idOf(typeTo1), _typeRegistry.get().idOf(typeTo2));
    auto top_pair = std::make_tuple(idOf(topologyTypeFrom1), idOf(topologyTypeFrom2));
    auto top_pair_to = std::make_tuple(idOf(topologyTypeTo1), idOf(topologyTypeTo2));
    addSpatialReaction(name, type_pair, top_pair, type_pair_to, top_pair_to, rate, radius, mode);
}

void TopologyRegistry::validateSpatialReaction(const spatial_reaction &reaction) const {
    if (reaction.rate() <= 0) {
        throw std::invalid_argument("The rate of an structural topology reaction (" + reaction.name() +
                                    ") should always be positive");
    }
    if (reaction.radius() <= 0) {
        throw std::invalid_argument("The radius of an structural topology reaction("
                                    + reaction.name() + ") should always be positive");
    }
    auto info1 = _typeRegistry.get().infoOf(reaction.type1());
    auto info2 = _typeRegistry.get().infoOf(reaction.type2());
    auto infoTo1 = _typeRegistry.get().infoOf(reaction.type_to1());
    auto infoTo2 = _typeRegistry.get().infoOf(reaction.type_to2());

    if (info1.flavor != particleflavor::TOPOLOGY && info2.flavor != particleflavor::TOPOLOGY) {
        throw std::invalid_argument(
                fmt::format("At least one of the educt types ({} and {}) of topology reaction {} need "
                                    "to be topology flavored!", info1.name, info2.name, reaction.name())
        );
    }

    {
        using tr_entry = spatial_reaction_map::value_type;
        auto it = std::find_if(_spatialReactions.begin(), _spatialReactions.end(),
                               [&reaction](const tr_entry &e) -> bool {
                                   return std::find_if(e.second.begin(), e.second.end(),
                                                       [&reaction](const spatial_reaction &tr) {
                                                           return tr.name() == reaction.name();
                                                       }) != e.second.end();
                               });
        if (it != _spatialReactions.end()) {
            throw std::invalid_argument("An structural topology reaction with the same name (" + reaction.name() +
                                        ") was already registered!");
        }
    }

    if (reaction.is_fusion()) {
        if (infoTo1.flavor != particleflavor::TOPOLOGY || infoTo2.flavor != particleflavor::TOPOLOGY) {
            throw std::invalid_argument(
                    fmt::format("One or both of the target types ({} and {}) of topology reaction {} did not have "
                                        "the topology flavor and is_fusion was set to true!",
                                infoTo1.name, infoTo2.name, reaction.name())
            );
        }
    } else {
        // flavors need to stay the same: topology particles must remain topology particles, same for normal particles
        if (info1.flavor != infoTo1.flavor) {
            throw std::invalid_argument(
                    fmt::format("if is_fusion is false, flavors need to remain same, which is not the case "
                                        "for {} (flavor={}) -> {} (flavor={})", info1.name,
                                particleflavor::particle_flavor_to_str(info1.flavor), infoTo1.name,
                                particleflavor::particle_flavor_to_str(infoTo1.flavor)));
        }
        if (info2.flavor != infoTo2.flavor) {
            throw std::invalid_argument(
                    fmt::format("if is_fusion is false, flavors need to remain same, which is not the case "
                                        "for {} (flavor={}) -> {} (flavor={})", info2.name,
                                particleflavor::particle_flavor_to_str(info2.flavor), infoTo2.name,
                                particleflavor::particle_flavor_to_str(infoTo2.flavor)));
        }
    }

}

void TopologyRegistry::configureBondPotential(const std::string &type1, const std::string &type2,
                                              const api::Bond &bond) {
    _potentialConfiguration.pairPotentials[std::make_tuple(_typeRegistry.get().idOf(type1),
                                                           _typeRegistry.get().idOf(type2))].push_back(bond);
}

void TopologyRegistry::configureAnglePotential(const std::string &type1, const std::string &type2,
                                               const std::string &type3, const api::Angle &angle) {
    _potentialConfiguration.anglePotentials[std::make_tuple(_typeRegistry.get().idOf(type1),
                                                            _typeRegistry.get().idOf(type2),
                                                            _typeRegistry.get().idOf(type3))].push_back(angle);
}

void TopologyRegistry::configureTorsionPotential(const std::string &type1, const std::string &type2,
                                                 const std::string &type3, const std::string &type4,
                                                 const api::TorsionAngle &torsionAngle) {
    _potentialConfiguration.torsionPotentials[std::make_tuple(_typeRegistry.get().idOf(type1),
                                                              _typeRegistry.get().idOf(type2),
                                                              _typeRegistry.get().idOf(type3),
                                                              _typeRegistry.get().idOf(type4))].push_back(
            torsionAngle);
}

}
}
}



