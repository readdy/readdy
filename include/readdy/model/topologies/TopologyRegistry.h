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
 * @file TopologyTypeRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 24.08.17
 * @copyright BSD-3
 */

#pragma once

#include <unordered_set>
#include <unordered_map>

#include <readdy/common/common.h>
#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/api/PotentialConfiguration.h>
#include <readdy/model/ParticleTypeRegistry.h>
#include <readdy/model/topologies/reactions/StructuralTopologyReaction.h>
#include <readdy/model/topologies/reactions/SpatialTopologyReaction.h>

#include "TopologyParticleTypeMap.h"

namespace readdy::model { class Context; }

namespace readdy::model::top {

struct TopologyType {
    using StructuralReaction = reactions::StructuralTopologyReaction;
    using StructuralReactionCollection = std::vector<StructuralReaction>;

    TopologyType() = default;

    TopologyType(std::string name, TopologyTypeId type, StructuralReactionCollection collection)
            : name(std::move(name)), type(type), structuralReactions(std::move(collection)) {}

    std::string name{""};
    TopologyTypeId type{0};
    StructuralReactionCollection structuralReactions{};
};

class TopologyRegistry {
public:

    using SpatialReaction = reactions::SpatialTopologyReaction;
    using SpatialReactionMap = topology_particle_type_tuple_umap<std::vector<SpatialReaction>>;
    using SpatialReactionCollection = SpatialReactionMap::mapped_type;
    using SpatialReactionTypeIds = std::unordered_set<ParticleTypeId>;

    using StructuralReaction = TopologyType::StructuralReaction;
    using StructuralReactionCollection = TopologyType::StructuralReactionCollection;

    using TypeCollection = std::vector<TopologyType>;

    TopologyTypeId addType(const std::string &name, const StructuralReactionCollection &reactions = {});

    void addStructuralReaction(TopologyTypeId typeId, reactions::StructuralTopologyReaction reaction) {
        {
            // test uniqueness of name
            for (const auto &topologyType : _registry) {
                for (const auto &structuralReaction : topologyType.structuralReactions) {
                    if (reaction.name() == structuralReaction.name()) {
                        throw std::invalid_argument(fmt::format("Structural topology reaction with name "
                                                                "\"{}\" already exists.", reaction.name()));
                    }
                }
            }
        }
        typeById(typeId).structuralReactions.push_back(std::move(reaction));
    }

    void addStructuralReaction(const std::string &type, const reactions::StructuralTopologyReaction &reaction) {
        addStructuralReaction(idOf(type), reaction);
    }

    void addStructuralReaction(const std::string &type, reactions::StructuralTopologyReaction &&reaction) {
        addStructuralReaction(idOf(type), std::forward<reactions::StructuralTopologyReaction>(reaction));
    }

    TopologyType &typeById(TopologyTypeId type) {
        auto it = std::find_if(std::begin(_registry), std::end(_registry), [type](const auto &e) {
            return e.type == type;
        });
        if (it != _registry.end()) {
            return *it;
        }
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }

    const TopologyType &typeById(TopologyTypeId type) const {
        auto it = std::find_if(std::begin(_registry), std::end(_registry), [type](const auto &e) {
            return e.type == type;
        });
        if (it != _registry.end()) {
            return *it;
        }
        throw std::invalid_argument(fmt::format("the requested type {} was not registered", type));
    }

    TopologyType &typeByName(const std::string &name) { return typeById(idOf(name)); }

    const TopologyType &typeByName(const std::string &name) const { return typeById(idOf(name)); }

    TypeCollection &types() { return _registry; }

    const TypeCollection &types() const { return _registry; }

    const std::string &nameOf(TopologyTypeId type) const {
        auto it = std::find_if(std::begin(_registry), std::end(_registry), [type](const auto &e) {
            return e.type == type;
        });
        if (it != _registry.end()) {
            return it->name;
        }
        throw std::invalid_argument(fmt::format("The requested type id {} did not exist.", type));
    }

    TopologyTypeId idOf(const std::string &name) const {
        if (name.empty()) return EmptyTopologyId;
        using entry_type = TypeCollection::value_type;
        auto it = std::find_if(std::begin(_registry), std::end(_registry), [&name](const auto &e) {
            return e.name == name;
        });
        if (it != _registry.end()) {
            return it->type;
        }
        throw std::invalid_argument(fmt::format("The requested type \"{}\" did not exist.", name));
    }

    reactions::SpatialTopologyReaction spatialTopologyReactionById(ReactionId id) const {
        for (const auto &[_, reactions] : _spatialReactions) {
            for (const auto &reaction : reactions) {
                if (reaction.id() == id) {
                    return reaction;
                }
            }
        }
        throw std::invalid_argument(fmt::format("Could not find spatial reaction with id {}.", id));
    }

    std::string spatialDescriptorById(ReactionId id) const {
        auto reaction = spatialTopologyReactionById(id);
        return generateSpatialReactionRepresentation(reaction);
    }

    std::string_view structuralNameById(ReactionId id) const {
        for (const auto &type : _registry) {
            for (const auto &structuralReaction : type.structuralReactions) {
                if (structuralReaction.id() == id) {
                    return structuralReaction.name();
                }
            }
        }
        throw std::invalid_argument(fmt::format("Could not find structural reaction with id {}.", id));
    }

    bool empty() {
        return _registry.empty();
    }

    const TopologyType::StructuralReactionCollection &structuralReactionsOf(TopologyTypeId type) const {
        auto it = std::find_if(_registry.begin(), _registry.end(), [type](const auto &e) {
            return e.type == type;
        });
        if (it != _registry.end()) {
            return it->structuralReactions;
        }
        throw std::invalid_argument(
                fmt::format("requested structural topology reactions of type {} which did not exist!", type)
        );
    }

    const StructuralReactionCollection &structuralReactionsOf(const std::string &type) const {
        return structuralReactionsOf(idOf(type));
    }

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
                               [](std::size_t current, const TypeCollection::const_iterator::value_type &other) {
                                   return current + other.structuralReactions.size();
                               }
        );
    }

    std::string generateSpatialReactionRepresentation(const SpatialReaction &reaction) const;

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

    void addSpatialReaction(const std::string &descriptor,
                            std::function<scalar(const GraphTopology &, const GraphTopology &)> rate_function,
                            scalar radius);

    void validateSpatialReaction(const SpatialReaction &reaction) const;

    const SpatialReactionMap &spatialReactionRegistry() const {
        return _spatialReactions;
    }

    const SpatialReactionCollection &spatialReactionsByType(
            const std::string &particleType1, const std::string &topologyType1,
            const std::string &particleType2, const std::string &topologyType2) const {
        return spatialReactionsByType(
                _types->idOf(particleType1), idOf(topologyType1),
                _types->idOf(particleType2), idOf(topologyType2));
    }

    const SpatialReactionCollection &spatialReactionsByType(ParticleTypeId t1, TopologyTypeId tt1,
                                                            ParticleTypeId t2, TopologyTypeId tt2) const {
        auto it = _spatialReactions.find(std::make_tuple(t1, tt1, t2, tt2));
        return it != _spatialReactions.end() ? it->second : _defaultTopologyReactions;
    }

    const SpatialReaction &spatialReactionByName(const std::string &name) const {
        for (const auto &e : _spatialReactions) {
            for (const auto &reaction : e.second) {
                if (reaction.name() == name) {
                    return reaction;
                }
            }
        }
        throw std::invalid_argument("No reaction with name \"" + name + "\" registered.");
    }

    SpatialReaction &spatialReactionByName(const std::string &name) {
        for (auto &e : _spatialReactions) {
            for (auto &reaction : e.second) {
                if (reaction.name() == name) {
                    return reaction;
                }
            }
        }
        throw std::invalid_argument("No reaction with name \"" + name + "\" registered.");
    }

    bool isSpatialReactionType(const std::string &name) const {
        return isSpatialReactionType(_types->idOf(name));
    }

    bool isSpatialReactionType(ParticleTypeId type) const {
        return _spatialReactionTypes.find(type) != _spatialReactionTypes.end();
    }

    /**
     * Configures a bond potential between a pair of particles with the given types. In order for the potential to come
     * into effect, the particles have to be connected in the connectivity graph of their topology. The types can be
     * reversed in order.
     * @param type1 first type of the tuple
     * @param type2 second type of the tuple
     * @param forceConstant the stiffness of the potential
     * @param length preferred distance between the particles with respect to this bond
     * @param type type of bond, defaults to harmonic bond
     */
    void configureBondPotential(const std::string &type1, const std::string &type2, const api::Bond &bond);

    /**
     * Configures an angle potential between a triple of particles with the given types. In order for the potential to
     * come into effect, the particles have to be connceted in the connectivity graph of their topology. The types
     * can be reversed in order.
     * @param type1 first type of the triple
     * @param type2 second type of the triple, type of the middle particle
     * @param type3 third type of the triple
     * @param forceConstant stiffness of the potential
     * @param equilibriumAngle the preferred angle between particles with respect to this potential
     * @param type the type of angle potential, defaults to harmonic angle
     */
    void configureAnglePotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                 const api::Angle &angle);

    /**
     * Configures a torsion potential between a quadruple of particles with the given types. In order for the potential
     * to come into effect, the particles have to be connected in the connectivity graph of their topology. The types
     * can be reversed in order.
     * @param type1 first type of the quadruple
     * @param type2 second type of the quadruple
     * @param type3 third type of the quadruple
     * @param type4 fourth type of the quadruple
     * @param forceConstant stiffness of the potential
     * @param multiplicity number of minima in the energy function
     * @param phi_0 reference torsion angle
     * @param type the type of torsion potential, defaults to cosine dihedral
     */
    void configureTorsionPotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                   const std::string &type4, const api::TorsionAngle &torsionAngle);

    api::PotentialConfiguration &potentialConfiguration() {
        return _potentialConfiguration;
    }

    const api::PotentialConfiguration &potentialConfiguration() const {
        return _potentialConfiguration;
    }

    const ParticleTypeRegistry &particleTypeRegistry() const {
        return *_types;
    }

private:
    static TopologyTypeId counter;

    TypeCollection _registry{};

    SpatialReactionMap _spatialReactions{};
    SpatialReactionCollection _defaultTopologyReactions{};
    SpatialReactionTypeIds _spatialReactionTypes{};

    const ParticleTypeRegistry *_types;

    api::PotentialConfiguration _potentialConfiguration{};

    friend readdy::model::Context;
};

}
