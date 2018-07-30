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
 * @file ReactionRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright BSD-3
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

    Reaction* order1ByName(const std::string &name) {
        for (auto &mapEntry : _o1Reactions) {
            for (auto &reaction : mapEntry.second) {
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

    Reaction* order2ByName(const std::string &name) {
        for (auto &mapEntry : _o2Reactions) {
            for (auto &reaction : mapEntry.second) {
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