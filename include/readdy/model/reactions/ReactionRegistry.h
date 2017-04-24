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
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class ReactionRegistry {
    using particle_t = readdy::model::Particle;
    using rea_ptr_vec1 = std::vector<std::unique_ptr<reactions::Reaction<1>>>;
    using rea_ptr_vec2 = std::vector<std::unique_ptr<reactions::Reaction<2>>>;
    using reaction_o1_registry_internal = std::unordered_map<particle_t::type_type, rea_ptr_vec1>;
    using reaction_o2_registry_internal = std::unordered_map<util::particle_type_pair, rea_ptr_vec2, util::particle_type_pair_hasher, util::particle_type_pair_equal_to>;

public:

    ReactionRegistry(std::reference_wrapper<const ParticleTypeRegistry> ref);

    ReactionRegistry(const ReactionRegistry &) = delete;

    ReactionRegistry &operator=(const ReactionRegistry &) = delete;

    ReactionRegistry(ReactionRegistry &&) = delete;

    ReactionRegistry &operator=(ReactionRegistry &&) = delete;

    using reaction_o1_registry = std::unordered_map<particle_t::type_type, std::vector<reactions::Reaction<1> *>>;
    using reaction_o2_registry = std::unordered_map<util::particle_type_pair, std::vector<reactions::Reaction<2> *>, util::particle_type_pair_hasher, util::particle_type_pair_equal_to>;

    const std::size_t &n_order1() const;

    const std::vector<const reactions::Reaction<1> *> order1_flat() const;

    const reactions::Reaction<1> *const order1_by_name(const std::string &name) const;

    const std::vector<reactions::Reaction<1> *> &order1_by_type(const particle_t::type_type type) const;

    const std::size_t &n_order2() const;

    const reaction_o2_registry &order2() const;

    const std::vector<const reactions::Reaction<2> *> order2_flat() const;

    const reactions::Reaction<2> *const order2_by_name(const std::string &name) const;

    const std::vector<reactions::Reaction<2> *> &order2_by_type(const particle_t::type_type type1,
                                                                const particle_t::type_type type2) const;


    const std::vector<reactions::Reaction<1> *> &order1_by_type(const std::string &type) const;

    const std::vector<reactions::Reaction<2> *> &order2_by_type(const std::string &type1,
                                                                const std::string &type2) const;

    template<typename R>
    const short add(std::unique_ptr<R> r,
                    typename std::enable_if<std::is_base_of<reactions::Reaction<1>, R>::value>::type * = 0) {
        log::trace("registering reaction {}", *r);
        const auto id = r->getId();
        const auto type = r->getEducts()[0];
        if (one_educt_registry_internal.find(type) == one_educt_registry_internal.end()) {
            one_educt_registry_internal.emplace(type, rea_ptr_vec1());
        }
        one_educt_registry_internal[type].push_back(std::move(r));
        n_order1_ += 1;
        return id;
    }

    template<typename R>
    const short add(std::unique_ptr<R> r,
                    typename std::enable_if<std::is_base_of<reactions::Reaction<2>, R>::value>::type * = 0) {
        log::trace("registering reaction {}", *r);
        const auto id = r->getId();
        const auto t1 = r->getEducts()[0];
        const auto t2 = r->getEducts()[1];

        const auto pp = std::tie(t1, t2);
        if (two_educts_registry_internal.find(pp) == two_educts_registry_internal.end()) {
            two_educts_registry_internal.emplace(pp, rea_ptr_vec2());
        }
        two_educts_registry_internal[pp].push_back(std::move(r));
        n_order2_ += 1;
        return id;
    }

    const short add_external(reactions::Reaction<1> *r);

    const short add_external(reactions::Reaction<2> *r);

    void configure();

    void debug_output() const;

private:
    using reaction_o1_registry_external = reaction_o1_registry;
    using reaction_o2_registry_external = reaction_o2_registry;

    std::size_t n_order1_ = 0;
    std::size_t n_order2_ = 0;

    const ParticleTypeRegistry &typeRegistry;

    reaction_o1_registry one_educt_registry{};
    reaction_o1_registry_internal one_educt_registry_internal{};
    reaction_o1_registry_external one_educt_registry_external{};
    reaction_o2_registry two_educts_registry{};
    reaction_o2_registry_internal two_educts_registry_internal{};
    reaction_o2_registry_external two_educts_registry_external{};

    std::vector<reactions::Reaction<1> *> defaultReactionsO1{};
    std::vector<reactions::Reaction<2> *> defaultReactionsO2{};
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)