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
public:
    using ReactionId = Reaction::ReactionId;

    using ReactionsCollection = std::vector<Reaction *>;
    using ReactionsO1Map = std::unordered_map<ParticleTypeId, ReactionsCollection>;
    using ReactionsO2Map = util::particle_type_pair_unordered_map<ReactionsCollection>;

    explicit ReactionRegistry(std::reference_wrapper<const ParticleTypeRegistry> ref) : _types(ref) {};

    ReactionRegistry(const ReactionRegistry &) = default;

    ReactionRegistry &operator=(const ReactionRegistry &) = default;

    ReactionRegistry(ReactionRegistry &&) = default;

    ReactionRegistry &operator=(ReactionRegistry &&) = default;

    ~ReactionRegistry() = default;

    const std::size_t &nOrder1() const {
        return _n_order1;
    }

    const ReactionsO1Map &order1() const {
        return _o1Reactions;
    }

    const ReactionsCollection order1Flat() const {
        ReactionsO1Map::mapped_type result;
        for (const auto &mapEntry : _o1Reactions) {
            for (const auto reaction : mapEntry.second) {
                result.push_back(reaction);
            }
        }
        return result;
    }

    const Reaction* order1ByName(const std::string &name) const {
        for (const auto &mapEntry : _o1Reactions) {
            for (const auto &reaction : mapEntry.second) {
                if (reaction->name() == name) return reaction;
            }
        }
        throw std::invalid_argument(fmt::format("No first order reaction with name \"{}\" found.", name));
    }

    const ReactionsCollection &order1ByType(const ParticleTypeId type) const {
        return readdy::util::collections::getOrDefault(_o1Reactions, type, DEFAULT_REACTIONS);
    }

    const std::size_t &nOrder2() const {
        return _n_order2;
    }

    const ReactionsO2Map &order2() const {
        return _o2Reactions;
    }

    const ReactionsCollection order2Flat() const {
        ReactionsO2Map::mapped_type result;
        for (const auto &mapEntry : _o2Reactions) {
            for (const auto reaction : mapEntry.second) {
                result.push_back(reaction);
            }
        }
        return result;
    }

    const Reaction* order2ByName(const std::string &name) const {
        for (const auto &mapEntry : _o2Reactions) {
            for (const auto &reaction : mapEntry.second) {
                if (reaction->name() == name) return reaction;
            }
        }
        throw std::invalid_argument(fmt::format("No second order reaction with name \"{}\" found.", name));
    }

    const ReactionsCollection &order2ByType(const ParticleTypeId type1, const ParticleTypeId type2) const {
        auto it = _o2Reactions.find(std::tie(type1, type2));
        return it != _o2Reactions.end() ? it->second : DEFAULT_REACTIONS;
    }

    const ReactionsCollection &order1ByType(const std::string &type) const {
        return order1ByType(_types(type));
    }

    const ReactionsCollection &order2ByType(const std::string &type1, const std::string &type2) const {
        return order2ByType(_types(type1), _types(type2));
    }

    std::string nameOf(ReactionId id) const;

    ReactionId idOf(const std::string &name) const;

    const Reaction* byId(ReactionId id) const;

    ReactionId add(const std::string &descriptor, scalar rate);

    /**
     * Method to register a conversion reaction "A->B".
     * @param name the name of the reaction
     * @param from the type of A
     * @param to the type of B
     * @param rate the rate at which this reaction is to be performed
     * @return a uuid
     */
    ReactionId addConversion(const std::string &name, const std::string &from, const std::string &to, scalar rate) {
        return addConversion(name, _types(from), _types(to), rate);
    }
    ReactionId addConversion(const std::string &name, ParticleTypeId from, ParticleTypeId to, scalar rate) {
        return emplaceReaction(std::make_shared<Conversion>(name, from, to, rate));
    }

    /**
     * Method to register an enzymatic reaction "A+C->B+C".
     * @param name the name of the reaction
     * @param catalyst the type of C
     * @param from the type of A
     * @param to the type of B
     * @param rate the rate at which this reaction is to be performed
     * @param eductDistance the distance at which B should be placed from C
     * @return a uuid
     */
    ReactionId addEnzymatic(const std::string &name, const std::string &catalyst, const std::string &from,
                             const std::string &to, scalar rate, scalar eductDistance) {
        return addEnzymatic(name, _types(catalyst), _types(from), _types(to),
                            rate, eductDistance);
    }
    ReactionId addEnzymatic(const std::string &name, ParticleTypeId catalyst, ParticleTypeId from,
                             ParticleTypeId to, scalar rate, scalar eductDistance) {
        return emplaceReaction(std::make_shared<Enzymatic>(name, catalyst, from, to, rate, eductDistance));
    }

    /**
     * Method to register a fission reaction "A->B+C".
     * @param name the name of the reaction
     * @param from the type of A
     * @param to1 the type of B
     * @param to2 the type of C
     * @param rate the rate at which this reaction is to be performed
     * @param productDistance the distance at which the products are placed
     * @param weight1 the weight for particle B with respect to the product distance
     * @param weight2 the weight for particle C with respect to the product distance
     * @return a uuid
     */
    ReactionId addFission(const std::string &name, const std::string &from, const std::string &to1,
                           const std::string &to2, scalar rate, scalar productDistance,
                           scalar weight1 = 0.5, scalar weight2 = 0.5) {
        return addFission(name, _types(from), _types(to1), _types(to2), rate,
                          productDistance, weight1, weight2);
    }
    ReactionId addFission(const std::string &name, ParticleTypeId from, ParticleTypeId to1,
                           ParticleTypeId to2, scalar rate, scalar productDistance,
                           scalar weight1 = 0.5, scalar weight2 = 0.5) {
        return emplaceReaction(std::make_shared<Fission>(name, from, to1, to2, rate, productDistance, weight1, weight2));
    }

    /**
     * Method to register a fusion reaction "A+B->C".
     * @param name the name of the reaction
     * @param from1 the type of A
     * @param from2 the type of B
     * @param to the type of C
     * @param rate the rate at which this reaction is to be performed
     * @param eductDistance the distance at which particles A and B become reactive
     * @param weight1 the weight of A with respect to the placement of C
     * @param weight2 the weight of B with respect to the placement of C
     * @return a uuid
     */
    ReactionId addFusion(const std::string &name, const std::string &from1, const std::string &from2,
                          const std::string &to, scalar rate, scalar eductDistance,
                          scalar weight1 = 0.5, scalar weight2 = 0.5) {
        return addFusion(name, _types(from1), _types(from2), _types(to), rate,
                         eductDistance, weight1, weight2);
    }
    ReactionId addFusion(const std::string &name, ParticleTypeId from1, ParticleTypeId from2,
                          ParticleTypeId to, scalar rate, scalar eductDistance,
                          scalar weight1 = 0.5, scalar weight2 = 0.5) {
        return emplaceReaction(std::make_shared<Fusion>(name, from1, from2, to, rate, eductDistance, weight1, weight2));
    }

    /**
     * Method to register a decay reaction.
     * @param name the name of the reaction
     * @param particleType the type for which this decay should be performed
     * @param rate the rate
     * @return a uuid
     */
    ReactionId addDecay(const std::string &name, const std::string &type, scalar rate) {
        return addDecay(name, _types(type), rate);
    }
    ReactionId addDecay(const std::string &name, ParticleTypeId type, scalar rate) {
        return emplaceReaction(std::make_shared<Decay>(name, type, rate));
    }

    std::string describe() const;

private:

    using OwnReactions = std::vector<std::shared_ptr<Reaction>>;
    using OwnReactionsO1Map = std::unordered_map<ParticleTypeId, OwnReactions>;
    using OwnReactionsO2Map = util::particle_type_pair_unordered_map<OwnReactions>;

    bool reactionNameExists(const std::string &name) const;

    ReactionRegistry::ReactionId emplaceReaction(const std::shared_ptr<Reaction> &reaction);

    std::size_t _n_order1{0};
    std::size_t _n_order2{0};

    std::reference_wrapper<const ParticleTypeRegistry> _types;

    ReactionsO1Map _o1Reactions{};
    ReactionsO2Map _o2Reactions{};

    OwnReactionsO1Map _ownO1Reactions{};
    OwnReactionsO2Map _ownO2Reactions{};

    static const ReactionsCollection DEFAULT_REACTIONS;
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)