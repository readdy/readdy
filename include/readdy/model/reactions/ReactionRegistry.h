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
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class ReactionRegistry {
    using particle = readdy::model::Particle;
    using rea_ptr_vec1 = std::vector<std::shared_ptr<reactions::Reaction<1>>>;
    using rea_ptr_vec2 = std::vector<std::shared_ptr<reactions::Reaction<2>>>;
    using reaction_o1_registry_internal = std::unordered_map<particle::type_type, rea_ptr_vec1>;
    using reaction_o2_registry_internal = util::particle_type_pair_unordered_map<rea_ptr_vec2>;

public:
    using reaction_o1_registry = std::unordered_map<particle::type_type, std::vector<reactions::Reaction<1> *>>;
    using reaction_o2_registry = util::particle_type_pair_unordered_map<std::vector<reactions::Reaction<2> *>>;
    using reaction_o2_types = std::unordered_set<particle_type_type>;

    using reactions_o1 = reaction_o1_registry::mapped_type;
    using reactions_o2 = reaction_o2_registry::mapped_type;
    using reaction_o1 = reactions_o1::value_type;
    using reaction_o2 = reactions_o2::value_type;

    using reaction_id = short;

    explicit ReactionRegistry(std::reference_wrapper<const ParticleTypeRegistry> ref);

    ReactionRegistry(const ReactionRegistry &) = default;

    ReactionRegistry &operator=(const ReactionRegistry &) = default;

    ReactionRegistry(ReactionRegistry &&) = default;

    ReactionRegistry &operator=(ReactionRegistry &&) = default;

    ~ReactionRegistry() = default;

    const std::size_t &n_order1() const;

    const reactions_o1 order1_flat() const;

    const reaction_o1 order1_by_name(const std::string &name) const;

    const reactions_o1 &order1_by_type(const particle::type_type type) const;

    const std::size_t &n_order2() const;

    const reaction_o2_registry &order2() const;

    const reactions_o2 order2_flat() const;

    const reaction_o2 order2_by_name(const std::string &name) const;

    const reactions_o2 &order2_by_type(const particle::type_type type1, const particle::type_type type2) const;

    const reactions_o1 &order1_by_type(const std::string &type) const;

    const reactions_o2 &order2_by_type(const std::string &type1, const std::string &type2) const;

    bool is_reaction_order2_type(particle_type_type type) const;

    reaction_id add(const std::string &descriptor, scalar rate);

    reaction_id addConversion(const std::string &name, const std::string &from, const std::string &to, scalar rate);

    reaction_id addConversion(const std::string &name, particle_type_type from, particle_type_type to, scalar rate);

    reaction_id addEnzymatic(const std::string &name, const std::string &catalyst, const std::string &from,
                             const std::string &to, scalar rate, scalar eductDistance);

    reaction_id addEnzymatic(const std::string &name, particle_type_type catalyst, particle_type_type from,
                             particle_type_type to, scalar rate, scalar eductDistance);

    reaction_id addFission(const std::string &name, const std::string &from, const std::string &to1,
                           const std::string &to2, scalar rate, scalar productDistance,
                           scalar weight1 = 0.5, scalar weight2 = 0.5);

    reaction_id addFission(const std::string &name, particle_type_type from, particle_type_type to1,
                           particle_type_type to2, scalar rate, scalar productDistance,
                           scalar weight1 = 0.5, scalar weight2 = 0.5);

    reaction_id addFusion(const std::string &name, const std::string &from1, const std::string &from2,
                          const std::string &to, scalar rate, scalar eductDistance,
                          scalar weight1 = 0.5, scalar weight2 = 0.5);

    reaction_id addFusion(const std::string &name, particle_type_type from1, particle_type_type from2,
                          particle_type_type to, scalar rate, scalar eductDistance,
                          scalar weight1 = 0.5, scalar weight2 = 0.5);

    reaction_id addDecay(const std::string &name, const std::string &type, scalar rate);

    reaction_id addDecay(const std::string &name, particle_type_type type, scalar rate);

    const short add_external(reaction_o1 r);

    const short add_external(reaction_o2 r);

    void configure();

    void debug_output() const;

private:
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

    reactions_o1 defaultReactionsO1{};
    reactions_o2 defaultReactionsO2{};
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)