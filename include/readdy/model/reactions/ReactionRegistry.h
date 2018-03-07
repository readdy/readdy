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
 * @file ReactionRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/model/ParticleTypeRegistry.h>
#include <unordered_set>
#include <readdy/common/Utils.h>
#include "Reaction.h"
#include "Enzymatic.h"
#include "Conversion.h"
#include "Fission.h"
#include "Fusion.h"
#include "Decay.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class ReactionRegistry {
    using particle = readdy::model::Particle;
    using reaction = Reaction;
    using rea_ptr_vec = std::vector<std::shared_ptr<reaction>>;
    using rea_raw_ptr_vec = std::vector<reaction *>;
    using reaction_o1_registry_internal = std::unordered_map<particle::type_type, rea_ptr_vec>;
    using reaction_o2_registry_internal = util::particle_type_pair_unordered_map<rea_ptr_vec>;

public:
    using reaction_o1_registry = std::unordered_map<particle::type_type, rea_raw_ptr_vec>;
    using reaction_o2_registry = util::particle_type_pair_unordered_map<rea_raw_ptr_vec>;
    using reaction_o2_types = std::unordered_set<particle_type_type>;

    using reactions = rea_raw_ptr_vec;

    using reaction_id = Reaction::reaction_id;
    using reactions_raw_ptr_map = std::unordered_map<reaction_id, reaction*>;

    explicit ReactionRegistry(std::reference_wrapper<const ParticleTypeRegistry> ref) : _types(ref) {};

    ReactionRegistry(const ReactionRegistry &) = default;

    ReactionRegistry &operator=(const ReactionRegistry &) = default;

    ReactionRegistry(ReactionRegistry &&) = default;

    ReactionRegistry &operator=(ReactionRegistry &&) = default;

    ~ReactionRegistry() = default;

    const std::size_t &nOrder1() const {
        return _n_order1;
    }

    const reaction_o1_registry &order1() const {
        return one_educt_registry;
    }

    const reactions order1Flat() const {
        reaction_o1_registry::mapped_type result;
        for (const auto &mapEntry : one_educt_registry) {
            for (const auto reaction : mapEntry.second) {
                result.push_back(reaction);
            }
        }
        return result;
    }

    const reaction* order1ByName(const std::string &name) const {
        for (const auto &mapEntry : one_educt_registry) {
            for (const auto &reaction : mapEntry.second) {
                if (reaction->name() == name) return reaction;
            }
        }
        throw std::invalid_argument(fmt::format("No first order reaction with name \"{}\" found.", name));
    }

    const reactions &order1ByType(const particle::type_type type) const {
        return readdy::util::collections::getOrDefault(one_educt_registry, type, defaultReactions);
    }

    const std::size_t &nOrder2() const {
        return _n_order2;
    }

    const reaction_o2_registry &order2() const {
        return two_educts_registry;
    }

    const reactions order2Flat() const {
        reaction_o2_registry::mapped_type result;
        for (const auto &mapEntry : two_educts_registry) {
            for (const auto reaction : mapEntry.second) {
                result.push_back(reaction);
            }
        }
        return result;
    }

    const reaction* order2ByName(const std::string &name) const {
        for (const auto &mapEntry : two_educts_registry) {
            for (const auto &reaction : mapEntry.second) {
                if (reaction->name() == name) return reaction;
            }
        }
        throw std::invalid_argument(fmt::format("No second order reaction with name \"{}\" found.", name));
    }

    const reactions &order2ByType(const particle::type_type type1, const particle::type_type type2) const {
        auto it = two_educts_registry.find(std::tie(type1, type2));
        return it != two_educts_registry.end() ? it->second : defaultReactions;
    }

    const reactions &order1ByType(const std::string &type) const {
        return order1ByType(_types.get().idOf(type));
    }

    const reactions &order2ByType(const std::string &type1, const std::string &type2) const {
        return order2ByType(_types.get().idOf(type1), _types.get().idOf(type2));
    }

    bool isReactionOrder2Type(particle_type_type type) const {
        return _reaction_o2_types.find(type) != _reaction_o2_types.end();
    }

    std::string nameOf(reaction_id id) const;

    reaction_id idOf(const std::string &name) const;

    const reaction* byId(reaction_id id) const;

    reaction_id add(const std::string &descriptor, scalar rate);

    reaction_id addConversion(const std::string &name, const std::string &from, const std::string &to, scalar rate) {
        return addConversion(name, _types.get().idOf(from), _types.get().idOf(to), rate);
    }

    reaction_id addConversion(const std::string &name, particle_type_type from, particle_type_type to, scalar rate) {
        return emplaceReaction(std::make_shared<Conversion>(name, from, to, rate));
    }

    reaction_id addEnzymatic(const std::string &name, const std::string &catalyst, const std::string &from,
                             const std::string &to, scalar rate, scalar eductDistance) {
        return addEnzymatic(name, _types.get().idOf(catalyst), _types.get().idOf(from), _types.get().idOf(to),
                            rate, eductDistance);
    }

    reaction_id addEnzymatic(const std::string &name, particle_type_type catalyst, particle_type_type from,
                             particle_type_type to, scalar rate, scalar eductDistance) {
        return emplaceReaction(std::make_shared<Enzymatic>(name, catalyst, from, to, rate, eductDistance));
    }

    reaction_id addFission(const std::string &name, const std::string &from, const std::string &to1,
                           const std::string &to2, scalar rate, scalar productDistance,
                           scalar weight1 = 0.5, scalar weight2 = 0.5) {
        return addFission(name, _types.get().idOf(from), _types.get().idOf(to1), _types.get().idOf(to2), rate,
                          productDistance, weight1, weight2);
    }

    reaction_id addFission(const std::string &name, particle_type_type from, particle_type_type to1,
                           particle_type_type to2, scalar rate, scalar productDistance,
                           scalar weight1 = 0.5, scalar weight2 = 0.5) {
        return emplaceReaction(std::make_shared<Fission>(name, from, to1, to2, rate, productDistance, weight1, weight2));
    }

    reaction_id addFusion(const std::string &name, const std::string &from1, const std::string &from2,
                          const std::string &to, scalar rate, scalar eductDistance,
                          scalar weight1 = 0.5, scalar weight2 = 0.5) {
        return addFusion(name, _types.get().idOf(from1), _types.get().idOf(from2), _types.get().idOf(to), rate,
                         eductDistance, weight1, weight2);
    }

    reaction_id addFusion(const std::string &name, particle_type_type from1, particle_type_type from2,
                          particle_type_type to, scalar rate, scalar eductDistance,
                          scalar weight1 = 0.5, scalar weight2 = 0.5) {
        return emplaceReaction(std::make_shared<Fusion>(name, from1, from2, to, rate, eductDistance, weight1, weight2));
    }

    reaction_id addDecay(const std::string &name, const std::string &type, scalar rate) {
        return addDecay(name, _types.get().idOf(type), rate);
    }

    reaction_id addDecay(const std::string &name, particle_type_type type, scalar rate) {
        return emplaceReaction(std::make_shared<Decay>(name, type, rate));
    }

    const short addExternal(reaction* r);

    void configure();

    std::string describe() const;

private:

    bool reactionNameExists(const std::string &name) const;

    ReactionRegistry::reaction_id emplaceReaction(const std::shared_ptr<Reaction> &reaction);

    using reaction_o1_registry_external = reaction_o1_registry;
    using reaction_o2_registry_external = reaction_o2_registry;

    std::size_t _n_order1{0};
    std::size_t _n_order2{0};

    std::reference_wrapper<const ParticleTypeRegistry> _types;

    reaction_o1_registry one_educt_registry{};
    reaction_o1_registry_internal one_educt_registry_internal{};
    reaction_o1_registry_external one_educt_registry_external{};
    reaction_o2_registry two_educts_registry{};
    reaction_o2_registry_internal two_educts_registry_internal{};
    reaction_o2_registry_external two_educts_registry_external{};
    reaction_o2_types _reaction_o2_types{};

    reactions defaultReactions{};
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)