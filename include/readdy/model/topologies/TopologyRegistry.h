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
 * @file TopologyTypeRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 24.08.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/model/topologies/reactions/StructuralTopologyReaction.h>
#include <unordered_set>
#include <unordered_map>
#include <readdy/model/topologies/reactions/SpatialTopologyReaction.h>
#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/model/ParticleTypeRegistry.h>
#include <readdy/api/PotentialConfiguration.h>
#include "TopologyParticleTypeMap.h"


NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

struct TopologyTypeInfo {
    using structural_reaction = reactions::StructuralTopologyReaction;
    using structural_reaction_vector = std::vector<structural_reaction>;

    std::string name{""};
    topology_type_type type{0};
    structural_reaction_vector structural_reactions{};
};

class TopologyRegistry {
public:

    using spatial_reaction = reactions::SpatialTopologyReaction;
    using spatial_reaction_map = topology_particle_type_tuple_umap<std::vector<spatial_reaction>>;
    using spatial_reactions = spatial_reaction_map::mapped_type;
    using spatial_reaction_types = std::unordered_set<particle_type_type>;

    using structural_reaction = TopologyTypeInfo::structural_reaction;
    using structural_reactions = TopologyTypeInfo::structural_reaction_vector;

    using type_registry = std::unordered_map<topology_type_type, TopologyTypeInfo>;

    explicit TopologyRegistry(const ParticleTypeRegistry &typeRegistry);

    TopologyRegistry(const TopologyRegistry &) = default;

    TopologyRegistry &operator=(const TopologyRegistry &) = default;

    TopologyRegistry(TopologyRegistry &&) = default;

    TopologyRegistry &operator=(TopologyRegistry &&) = default;

    ~TopologyRegistry() = default;

    topology_type_type addType(const std::string &name, const structural_reactions &reactions = {});

    void addStructuralReaction(topology_type_type type, const reactions::StructuralTopologyReaction &reaction);

    void addStructuralReaction(topology_type_type type, reactions::StructuralTopologyReaction &&reaction);

    void addStructuralReaction(const std::string &type, const reactions::StructuralTopologyReaction &reaction) {
        addStructuralReaction(idOf(type), reaction);
    }

    void addStructuralReaction(const std::string &type, reactions::StructuralTopologyReaction &&reaction) {
        addStructuralReaction(idOf(type), std::forward<reactions::StructuralTopologyReaction>(reaction));
    }

    const TopologyTypeInfo &infoOf(topology_type_type type) const {
        auto it = _registry.find(type);
        if (it != _registry.end()) {
            return it->second;
        }
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }

    const std::string &nameOf(topology_type_type type) const {
        auto it = _registry.find(type);
        if (it != _registry.end()) {
            return it->second.name;
        }
        throw std::invalid_argument(fmt::format("The requested type id {} did not exist.", type));
    }

    topology_type_type idOf(const std::string &name) const {
        if (name.empty()) return topology_type_empty;
        using entry_type = type_registry::value_type;
        auto it = std::find_if(_registry.begin(), _registry.end(), [&name](const entry_type &entry) {
            return entry.second.name == name;
        });
        if (it != _registry.end()) {
            return it->first;
        }
        throw std::invalid_argument(fmt::format("The requested type \"{}\" did not exist.", name));
    }

    bool empty() {
        return _registry.empty();
    }

    bool containsStructuralReactions() {
        return _containsStructuralReactions;
    }

    const TopologyTypeInfo::structural_reaction_vector &structuralReactionsOf(topology_type_type type) const {
        auto it = _registry.find(type);
        if (it != _registry.end()) {
            return it->second.structural_reactions;
        }
        throw std::invalid_argument(fmt::format("requested structural topology reactions of type {} which did not exist!", type));
        log::warn("requested structural topology reactions of type {} which did not exist!", type);
        return _defaultInfo.structural_reactions;
    }

    const structural_reactions &structuralReactionsOf(const std::string &type) const {
        return structuralReactionsOf(idOf(type));
    }

    void configure();

    std::string describe() const;

    void addSpatialReaction(reactions::SpatialTopologyReaction &&reaction);

    void addSpatialReaction(const std::string &name, const std::string &typeFrom1,
                            const std::string &typeFrom2, const std::string &topologyTypeFrom1,
                            const std::string &topologyTypeFrom2, const std::string &typeTo1,
                            const std::string &typeTo2, const std::string &topologyTypeTo1,
                            const std::string &topologyTypeTo2, scalar rate, scalar radius, reactions::STRMode mode);

    void addSpatialReaction(const std::string &name, const util::particle_type_pair &types,
                            const topology_type_pair &topology_types,
                            const util::particle_type_pair &types_to, const topology_type_pair &topology_types_to,
                            scalar rate, scalar radius, reactions::STRMode mode);

    std::size_t nStructuralReactions() const {
        return std::accumulate(_registry.begin(), _registry.end(), 0_z,
                               [](std::size_t current, const type_registry::const_iterator::value_type &other) {
                                   return current + other.second.structural_reactions.size();
                               }
        );
    }

    std::string generateSpatialReactionRepresentation(const spatial_reaction &reaction) const;

    /**
     * Convenience function for adding spatial reaction.
     *
     * Pass descriptor of form for Topology<->Topology
     * - fusion type:
     *     name: T1 (p1) + T2 (p2) -> T3 (p3 + p4) [self=true|false]
     * - enzymatic type:
     *     name: T1 (p1) + T2 (p2) -> T3 (p3) + T4 (p4)
     *
     * and for Topology<->Particle
     * - fusion type:
     *     name: T1 (p1) + p2 -> T2 (p3 + p4)
     * - enzymatic type:
     *     name: T1 (p1) + p2 -> T2 (p3) + p4
     *
     * @param descriptor descriptor as described above
     * @param rate the rate
     * @param radius the radius
     */
    void addSpatialReaction(const std::string &descriptor, scalar rate, scalar radius);

    void validateSpatialReaction(const spatial_reaction &reaction) const;

    const spatial_reaction_map &spatialReactionRegistry() const {
        return _spatialReactions;
    }

    const spatial_reactions &spatialReactionsByType(particle_type_type t1, topology_type_type tt1,
                                                    particle_type_type t2, topology_type_type tt2) const {
        auto it = _spatialReactions.find(std::make_tuple(t1, tt1, t2, tt2));
        return it != _spatialReactions.end() ? it->second : _defaultTopologyReactions;
    }

    bool isSpatialReactionType(const std::string &name) const {
        return isSpatialReactionType(_typeRegistry.get().idOf(name));
    }

    bool isSpatialReactionType(particle_type_type type) const {
        return _topologyReactionTypes.find(type) != _topologyReactionTypes.end();
    }

    /*
     * Potentials
     */

    void configureBondPotential(const std::string &type1, const std::string &type2, const api::Bond &bond);

    void configureAnglePotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                 const api::Angle &angle);

    void configureTorsionPotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                   const std::string &type4, const api::TorsionAngle &torsionAngle);

    api::PotentialConfiguration &potentialConfiguration() {
        return _potentialConfiguration;
    }

    const api::PotentialConfiguration &potentialConfiguration() const {
        return _potentialConfiguration;
    }

    const ParticleTypeRegistry &particleTypeRegistry() const {
        return _typeRegistry.get();
    }

private:
    static topology_type_type counter;

    TopologyTypeInfo _defaultInfo;
    bool _containsStructuralReactions{false};
    type_registry _registry{};

    spatial_reaction_map _spatialReactions{};
    spatial_reactions _defaultTopologyReactions{};
    spatial_reaction_types _topologyReactionTypes{};

    std::reference_wrapper<const ParticleTypeRegistry> _typeRegistry;

    api::PotentialConfiguration _potentialConfiguration{};
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)