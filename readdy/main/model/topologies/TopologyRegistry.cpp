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
 * << detailed description >>
 *
 * @file #.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 24.08.17
 * @copyright GPL-3
 */

#include <readdy/model/topologies/TopologyRegistry.h>
#include <readdy/model/Utils.h>

namespace readdy {
namespace model {
namespace top {
TopologyTypeId TopologyRegistry::counter = 0;

void TopologyRegistry::addStructuralReaction(TopologyTypeId id,
                                             const reactions::StructuralTopologyReaction &reaction) {
    typeById(id).structuralReactions.push_back(reaction);
}

void TopologyRegistry::addStructuralReaction(TopologyTypeId id,
                                             reactions::StructuralTopologyReaction &&reaction) {
    typeById(id).structuralReactions.push_back(std::move(reaction));
}

std::string TopologyRegistry::describe() const {
    namespace rus = readdy::util::str;
    std::string description;

    if(!_potentialConfiguration.pairPotentials.empty()
       || !_potentialConfiguration.anglePotentials.empty()
       || !_potentialConfiguration.torsionPotentials.empty()) {
        description += fmt::format(" - topology potential configuration:{}", rus::newline);

        if(!_potentialConfiguration.pairPotentials.empty()) {
            description += fmt::format("     - bonds ({}):{}", _potentialConfiguration.pairPotentials.size(),
                                       rus::newline);
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
                    description += fmt::format("             * {} bond with force constant {} and length {}{}",
                                               bondToStr(bond),
                                               bond.forceConstant, bond.length, rus::newline);
                }
            }
        }

        if(!_potentialConfiguration.anglePotentials.empty()) {
            description += fmt::format("     - angles ({}):{}", _potentialConfiguration.anglePotentials.size(),
                                       rus::newline);
            for (const auto &entry : _potentialConfiguration.anglePotentials) {
                auto angleToStr = [](const api::Angle &angle) -> std::string {
                    switch (angle.type) {
                        case api::AngleType::HARMONIC:
                            return "Harmonic";
                    }
                };
                for (const auto &angle : entry.second) {
                    description += fmt::format(
                            "             * {} angle with force constant {} and equilibrium angle {}{}",
                            angleToStr(angle), angle.forceConstant, angle.equilibriumAngle,
                            rus::newline);
                }
            }
        }
        if(!_potentialConfiguration.torsionPotentials.empty()) {
            description += fmt::format("     - torsions ({}):{}", _potentialConfiguration.torsionPotentials.size(),
                                       rus::newline);
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
        }
    }

    if(!_registry.empty()) {
        description += fmt::format(" - topology types:{}", rus::newline);
        for (const auto &entry : _registry) {
            description += fmt::format("     * topology type \"{}\" with {} structural reactions{}",
                                       entry.name, entry.structuralReactions.size(), rus::newline);
        }
    }

    if (!_spatialReactions.empty()) {
        description += fmt::format(" - spatial topology reactions:{}", rus::newline);
        for (const auto &entry : _spatialReactions) {
            for (const auto &reaction : entry.second) {
                std::stringstream ss;
                switch(reaction.mode()) {
                    case reactions::STRMode::TT_ENZYMATIC: {
                        ss << "Topology-topology enzymatic reaction \"";
                        break;
                    }
                    case reactions::STRMode::TT_FUSION:
                    case reactions::STRMode::TT_FUSION_ALLOW_SELF: {
                        ss << "Topology-topology fusion reaction \"";
                        break;
                    }
                    case reactions::STRMode::TP_ENZYMATIC: {
                        ss << "Topology-particle enzymatic reaction \"";
                        break;
                    }
                    case reactions::STRMode::TP_FUSION: {
                        ss << "Topology-particle fusion reaction \"";
                        break;
                    }
                }
                ss << generateSpatialReactionRepresentation(reaction);
                ss << "\"";

                description += fmt::format("     * {}{}", ss.str(), rus::newline);
            }
        }
    }

    if(nStructuralReactions() > 0) {
        description += fmt::format(" - structural topology reactions:{}", rus::newline);
        for (const auto &entry : _registry) {
            if(!entry.structuralReactions.empty()) {
                description += fmt::format("     - for topology type \"{}\" with {} structural reactions:{}",
                                           entry.name, entry.structuralReactions.size(), rus::newline);
                for (const auto &r : entry.structuralReactions) {
                    description += fmt::format("         * reaction with roll_back = {} and create child topologies = {}{}",
                                               r.rolls_back_if_invalid(), r.creates_child_topologies_after_reaction(),
                                               rus::newline);
                }
            }
        }
    }
    return description;
}

readdy::TopologyTypeId TopologyRegistry::addType(const std::string &name,
                                                   const StructuralReactionCollection &reactions) {
    util::validateTypeName(name);
    using entry_type = TypeCollection::value_type;
    auto it = std::find_if(_registry.begin(), _registry.end(), [&name](const entry_type &entry) {
        return entry.name == name;
    });
    if (it == _registry.end()) {
        const auto type = counter++;
        _registry.emplace_back(name, type, reactions);
        return type;
    }
    throw std::invalid_argument(fmt::format("A type with name \"{}\" already existed.", name));
}

TopologyRegistry::TopologyRegistry(
        const readdy::model::ParticleTypeRegistry &typeRegistry) : _typeRegistry(typeRegistry) {}


void TopologyRegistry::addSpatialReaction(const std::string &name, const readdy::util::particle_type_pair &types,
                                          const topology_type_pair &topology_types,
                                          const readdy::util::particle_type_pair &types_to,
                                          const topology_type_pair &topology_types_to,
                                          scalar rate, scalar radius, reactions::STRMode mode) {
    SpatialReaction reaction(name, types, topology_types, types_to, topology_types_to, rate, radius, mode);
    addSpatialReaction(std::move(reaction));
}

void TopologyRegistry::addSpatialReaction(const std::string &descriptor, scalar rate, scalar radius) {
    reactions::STRParser parser(*this);
    addSpatialReaction(parser.parse(descriptor, rate, radius));
}

void TopologyRegistry::addSpatialReaction(reactions::SpatialTopologyReaction &&reaction) {
    validateSpatialReaction(reaction);
    auto key = std::make_tuple(reaction.type1(), reaction.top_type1(), reaction.type2(), reaction.top_type2());
    _spatialReactionTypes.emplace(reaction.type1());
    _spatialReactionTypes.emplace(reaction.type2());
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

void TopologyRegistry::validateSpatialReaction(const SpatialReaction &reaction) const {
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
        using tr_entry = SpatialReactionMap::value_type;
        auto it = std::find_if(_spatialReactions.begin(), _spatialReactions.end(),
                               [&reaction](const tr_entry &e) -> bool {
                                   return std::find_if(e.second.begin(), e.second.end(),
                                                       [&reaction](const SpatialReaction &tr) {
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

std::string TopologyRegistry::generateSpatialReactionRepresentation(const SpatialReaction &reaction) const {
    auto pName = [&](ParticleTypeId t) { return _typeRegistry.get().nameOf(t); };
    auto tName = [&](TopologyTypeId t) { return nameOf(t); };

    std::stringstream ss;
    ss << reaction.name() << ": ";
    switch (reaction.mode()) {
        case reactions::STRMode::TT_ENZYMATIC: {
            ss << fmt::format("{}({}) + {}({}) -> {}({}) + {}({})",
                              tName(reaction.top_type1()), pName(reaction.type1()),
                              tName(reaction.top_type2()), pName(reaction.type2()),
                              tName(reaction.top_type_to1()), pName(reaction.type_to1()),
                              tName(reaction.top_type_to2()), pName(reaction.type_to2()));
            break;
        }
        case reactions::STRMode::TT_FUSION:
        case reactions::STRMode::TT_FUSION_ALLOW_SELF: {
            ss << fmt::format("{}({}) + {}({}) -> {}({}--{})",
                              tName(reaction.top_type1()), pName(reaction.type1()),
                              tName(reaction.top_type2()), pName(reaction.type2()),
                              tName(reaction.top_type_to1()), pName(reaction.type_to1()), pName(reaction.type_to2()));
            if(reaction.allow_self_connection()) {
                ss << " [self=true]";
            }
            break;
        }
        case reactions::STRMode::TP_ENZYMATIC: {
            ss << fmt::format("{}({}) + ({}) -> {}({}) + ({})", tName(reaction.top_type1()), pName(reaction.type1()),
                              pName(reaction.type2()), tName(reaction.top_type_to1()), pName(reaction.type_to1()),
                              pName(reaction.type_to2()));
            break;
        }
        case reactions::STRMode::TP_FUSION: {
            ss << fmt::format("{}({}) + ({}) -> {}({}--{})", tName(reaction.top_type1()), pName(reaction.type1()),
                              pName(reaction.type2()), tName(reaction.top_type_to1()), pName(reaction.type_to1()),
                              pName(reaction.type_to2()));
            break;
        }
    }
    return ss.str();
}

}
}
}
