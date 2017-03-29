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
 * @file ReactionRegistry.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/reactions/ReactionRegistry.h>
#include <readdy/common/Utils.h>

namespace readdy {
namespace model {
namespace reactions {

const std::vector<const reactions::Reaction<1> *> ReactionRegistry::order1_flat() const {
    auto result = std::vector<const reactions::Reaction<1> *>();
    for (const auto &mapEntry : reactionOneEductRegistry) {
        for (const auto reaction : mapEntry.second) {
            result.push_back(reaction);
        }
    }
    return result;
}

const reactions::Reaction<1> *const ReactionRegistry::order1_by_name(const std::string &name) const {
    for (const auto &mapEntry : reactionOneEductRegistry) {
        for (const auto &reaction : mapEntry.second) {
            if (reaction->getName() == name) return reaction;
        }
    }

    return nullptr;
}

const std::vector<const reactions::Reaction<2> *> ReactionRegistry::order2_flat() const {
    auto result = std::vector<const reactions::Reaction<2> *>();
    for (const auto &mapEntry : reactionTwoEductsRegistry) {
        for (const auto reaction : mapEntry.second) {
            result.push_back(reaction);
        }
    }
    return result;
}

const std::vector<Reaction<1> *> &ReactionRegistry::order1_by_type(const Particle::type_type type) const {
    return util::collections::getOrDefault(reactionOneEductRegistry, type, defaultReactionsO1);
}

const reactions::Reaction<2> *const ReactionRegistry::order2_by_name(const std::string &name) const {
    for (const auto &mapEntry : reactionTwoEductsRegistry) {
        for (const auto &reaction : mapEntry.second) {
            if (reaction->getName() == name) return reaction;
        }
    }
    return nullptr;
}

const std::vector<Reaction<2> *> &
ReactionRegistry::order2_by_type(const Particle::type_type type1, const Particle::type_type type2) const {
    decltype(reactionTwoEductsRegistry.find(std::tie(type1, type2))) it;
    if((it = reactionTwoEductsRegistry.find(std::tie(type1, type2))) != reactionTwoEductsRegistry.end()) {
        return it->second;
    }
    return defaultReactionsO2;
}

void ReactionRegistry::configure() {
    namespace coll = readdy::util::collections;
    using pair = util::particle_type_pair;
    using reaction1ptr = std::unique_ptr<reactions::Reaction<1>>;
    using reaction2ptr = std::unique_ptr<reactions::Reaction<2>>;

    reactionOneEductRegistry.clear();
    reactionTwoEductsRegistry.clear();
    coll::for_each_value(reactionOneEductRegistryInternal,
                         [&](const particle_t::type_type type, const reaction1ptr &ptr) {
                             (reactionOneEductRegistry)[type].push_back(ptr.get());
                         });
    coll::for_each_value(reactionTwoEductsRegistryInternal, [&](const pair &type, const reaction2ptr &r) {
        (reactionTwoEductsRegistry)[type].push_back(r.get());
    });
    coll::for_each_value(reactionOneEductRegistryExternal,
                         [&](const particle_t::type_type type, reactions::Reaction<1> *ptr) {
                             (reactionOneEductRegistry)[type].push_back(ptr);
                         });
    coll::for_each_value(reactionTwoEductsRegistryExternal, [&](const pair &type, reactions::Reaction<2> *r) {
        (reactionTwoEductsRegistry)[type].push_back(r);
    });
}

void ReactionRegistry::debug_output() const {
    if (!reactionOneEductRegistry.empty()) {
        log::debug(" - reactions of order 1:");
        for(const auto& entry : reactionOneEductRegistry) {
            for(const auto& reaction : entry.second) {
                log::debug("     * reaction {}", *reaction);

            }
        }
    }
    if (!reactionTwoEductsRegistry.empty()) {
        log::debug(" - reactions of order 2:");
        for(const auto& entry : reactionTwoEductsRegistry) {
            for(const auto& reaction : entry.second) {
                log::debug("     * reaction {}", *reaction);
            }
        }
    }
}

const std::size_t &ReactionRegistry::n_order1() const {
    return n_order1_;
}

const std::size_t &ReactionRegistry::n_order2() const {
    return n_order2_;
}

const short ReactionRegistry::add_external(reactions::Reaction<1> *r) {
    reactionOneEductRegistryExternal[r->getEducts()[0]].push_back(r);
    ++n_order1_;
    return r->getId();
}

const short ReactionRegistry::add_external(reactions::Reaction<2> *r) {
    reactionTwoEductsRegistryExternal[std::tie(r->getEducts()[0], r->getEducts()[1])].push_back(r);
    ++n_order2_;
    return r->getId();
}

}
}
}