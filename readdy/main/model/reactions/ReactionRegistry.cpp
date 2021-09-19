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
 * @file ReactionRegistry.cpp
 * @brief << brief description >>
 * @author clonker
 * @author chrisfroe
 * @date 29.03.17
 * @copyright BSD-3
 */

#include <readdy/model/reactions/ReactionRegistry.h>
#include <readdy/common/Utils.h>
#include <readdy/model/reactions/Conversion.h>
#include <readdy/model/reactions/Enzymatic.h>
#include <readdy/model/reactions/Fission.h>
#include <readdy/model/reactions/Fusion.h>
#include <readdy/model/reactions/Decay.h>
#include <readdy/common/string.h>
#include <readdy/model/Utils.h>
#include <regex>
#include <utility>

namespace readdy::model::reactions {

ReactionId ReactionRegistry::emplaceReaction(const std::shared_ptr<Reaction> &reaction) {
    if (reactionNameExists(reaction->name())) {
        throw std::invalid_argument(fmt::format("A reaction with the name {} exists already", reaction->name()));
    }
    ReactionId id = reaction->id();
    if (reaction->nEducts() == 1) {
        auto t = reaction->educts()[0];
        _ownO1Reactions[t].push_back(reaction);
        _o1Reactions[t].push_back(_ownO1Reactions[t].back().get());
        _n_order1 += 1;
    } else {
        auto pp = std::make_tuple(reaction->educts()[0], reaction->educts()[1]);
        _ownO2Reactions[pp].push_back(reaction);
        _o2Reactions[pp].push_back(_ownO2Reactions[pp].back().get());
        _n_order2 += 1;
    }
    return id;
}

std::string ReactionRegistry::describe() const {
    namespace rus = readdy::util::str;
    std::string description;
    auto nameOf = [&](ParticleTypeId t) {return _types->nameOf(t);};

    if (!_o1Reactions.empty()) {
        description += fmt::format(" - unimolecular reactions:\n");
        for (const auto &entry : _o1Reactions) {
            for (const auto reaction : entry.second) {
                const auto &educts = reaction->educts();
                const auto &products = reaction->products();

                switch (reaction->type()) {
                    case ReactionType::Conversion: {
                        description += fmt::format("     * Conversion \"{}\": {} -> {} with a rate of {}\n",
                                                   reaction->name(), nameOf(educts[0]), nameOf(products[0]),
                                                   reaction->rate());
                        break;
                    }
                    case ReactionType::Fission: {
                        description += fmt::format("     * Fission \"{}\": {} -> {} + {} with a rate of {}, a product distance of {}, and weights {} and {}\n",
                                                   reaction->name(), nameOf(educts[0]), nameOf(products[0]),
                                                   nameOf(products[1]), reaction->rate(), reaction->productDistance(),
                                                   reaction->weight1(), reaction->weight2());
                        break;
                    }
                    case ReactionType::Decay: {
                        description += fmt::format("     * Decay \"{}\": {} -> ø with a rate of {}\n",
                                                   reaction->name(), nameOf(educts[0]), reaction->rate());
                        break;
                    }
                    default: {
                        throw std::logic_error(fmt::format("unknown unimolecular reaction in registry that "
                                                           "was of type {}", reaction->type()));
                    }
                }
            }
        }
    }
    if (!_o2Reactions.empty()) {
        description += fmt::format(" - bimolecular reactions:\n");
        for (const auto &entry : _o2Reactions) {
            for (const auto reaction : entry.second) {
                const auto &educts = reaction->educts();
                const auto &products = reaction->products();

                switch(reaction->type()) {
                    case ReactionType::Enzymatic: {
                        auto enzymatic = dynamic_cast<const Enzymatic*>(reaction);
                        description += fmt::format("     * Enzymatic \"{}\": {} + {} -> {} + {} with a rate of {} and an educt distance of {}\n",
                                                   enzymatic->name(), nameOf(enzymatic->getFrom()), nameOf(enzymatic->getCatalyst()),
                                                   nameOf(enzymatic->getTo()), nameOf(enzymatic->getCatalyst()),
                                                   enzymatic->rate(), enzymatic->eductDistance());
                        break;
                    }
                    case ReactionType::Fusion: {
                        auto fusion = dynamic_cast<const Fusion*>(reaction);
                        description += fmt::format("     * Fusion \"{}\": {} + {} -> {} with a rate of {}, an educt distance of {}, and weights {} and {}\n",
                                                   fusion->name(), nameOf(fusion->getFrom1()), nameOf(fusion->getFrom2()),
                                                   nameOf(fusion->getTo()), fusion->rate(), fusion->eductDistance(),
                                                   fusion->weight1(), fusion->weight2());
                        break;
                    }
                    default: {
                        throw std::logic_error(fmt::format("unknown bimolecular reaction in registry that "
                                                           "was of type {}", reaction->type()));
                    }
                }
            }
        }
    }
    return description;
}

ReactionId ReactionRegistry::add(const std::string &descriptor, scalar rate) {
    namespace mutil = readdy::model::util;
    namespace rutil = readdy::util;
    log::trace("begin parsing \"{}\"", descriptor);
    auto arrowPos = descriptor.find(mutil::arrow);
    if (arrowPos == std::string::npos) {
        throw std::invalid_argument(fmt::format(
                "the descriptor must contain an arrow (\"{}\") to indicate lhs and rhs.", mutil::arrow
        ));
    }
    if (descriptor.find(mutil::arrow, arrowPos + 1) != std::string::npos) {
        throw std::invalid_argument(fmt::format(
                "the descriptor must not contain more than one arrow (\"{}\").", mutil::arrow
        ));
    }
    auto lhs = descriptor.substr(0, arrowPos);
    auto rhs = descriptor.substr(arrowPos + std::strlen(mutil::arrow), std::string::npos);

    rutil::str::trim(lhs);
    rutil::str::trim(rhs);

    std::string name;
    {
        auto colonPos = lhs.find(':');
        if (colonPos == std::string::npos) {
            throw std::invalid_argument("The descriptor did not contain a colon ':' to specify the end of the name.");
        }
        name = rutil::str::trim_copy(lhs.substr(0, colonPos));
        lhs = rutil::str::trim_copy(lhs.substr(colonPos + 1, std::string::npos));
    }

    // some helper functions
    static std::regex parenthesesContentRegex(R"(\(([^\)]*)\))");
    static auto getReactionRadius = [](const std::string &s) {
        std::smatch radiusMatch;
        if (std::regex_search(s, radiusMatch, parenthesesContentRegex)) {
            auto radiusString = rutil::str::trim_copy(radiusMatch.str());
            return rutil::str::trim_copy(radiusString.substr(1, radiusString.size() - 2)); // remove parentheses
        }
        throw std::invalid_argument(fmt::format("The term \"{}\" did not contain a pair of parentheses '(*)'", s));
    };

    static std::regex bracketsContentRegex(R"(\[([^\)]*)\])");
    static auto getWeights = [](const std::string &s) {
        std::smatch weightsMatch;
        if (std::regex_search(s, weightsMatch, bracketsContentRegex)) {
            auto weightsString = rutil::str::trim_copy(weightsMatch.str());
            auto commaPos = weightsString.find(',');
            if (commaPos == std::string::npos) {
                throw std::invalid_argument(fmt::format("The term \"{}\" did not contain a comma ','", s));
            }
            auto w1 = rutil::str::trim_copy(weightsString.substr(1, commaPos));
            auto w2 = rutil::str::trim_copy(weightsString.substr(commaPos + 1, weightsString.size() - 2));
            return std::make_tuple(w1, w2);
        }
        throw std::invalid_argument(fmt::format("The term \"{}\" did not contain a pair of brackets '[*]'", s));
    };

    static auto treatSideWithPlus = [](const std::string &s, std::size_t plusPos, bool searchRadius = false) {
        auto pType1 = rutil::str::trim_copy(s.substr(0, plusPos));
        std::string pType2;
        std::string reactionRadius;
        if (searchRadius) {
            auto closingParentheses = s.find(')');
            if (closingParentheses == std::string::npos) {
                throw std::invalid_argument(fmt::format("The side (\"{}\") did not contain a ')'.", s));
            }
            pType2 = rutil::str::trim_copy(s.substr(closingParentheses + 1, std::string::npos));
            reactionRadius = getReactionRadius(s);
        } else {
            pType2 = rutil::str::trim_copy(s.substr(plusPos + 1, std::string::npos));
        }
        return std::make_tuple(pType1, pType2, reactionRadius);
    };

    std::string w1, w2;
    scalar weight1{0.5}, weight2{0.5};
    {
        auto bracketPos = rhs.find('[');
        if (bracketPos != std::string::npos) {
            std::tie(w1, w2) = getWeights(rhs);
            weight1 = static_cast<scalar>(std::stod(w1));
            weight2 = static_cast<scalar>(std::stod(w2));
            rhs = rutil::str::trim_copy(rhs.substr(0, bracketPos));
        }
    }

    std::size_t numberEducts;
    std::string educt1, educt2, eductDistance;
    {
        auto plusPos = lhs.find('+');
        if (plusPos == std::string::npos) {
            numberEducts = 1;
            educt1 = lhs;
        } else {
            numberEducts = 2;
            std::tie(educt1, educt2, eductDistance) = treatSideWithPlus(lhs, plusPos, true);
        }
    }

    std::string product1, product2, productDistance;
    {
        if (rhs.empty()) {
            // numberProducts = 0 -> decay
            if (numberEducts != 1) {
                throw std::invalid_argument(fmt::format("Left hand side (\"{}\") did not contain one educt", lhs));
            }
            return addDecay(name, educt1, rate);
        } else {
            auto plusPos = rhs.find('+');
            if (plusPos == std::string::npos) {
                // numberProducts = 1 -> either conversion or fusion
                product1 = rhs;
                if (numberEducts == 1) {
                    // conversion
                    return addConversion(name, educt1, product1, rate);
                } else {
                    // fusion
                    auto distance = static_cast<scalar>(std::stod(eductDistance));
                    return addFusion(name, educt1, educt2, product1, rate, distance, weight1, weight2);
                }
            } else {
                // numberProducts = 2 -> either fission or enzymatic
                if (numberEducts == 1) {
                    // fission
                    std::tie(product1, product2, productDistance) = treatSideWithPlus(rhs, plusPos, true);
                    auto distance = static_cast<scalar>(std::stod(productDistance));
                    return addFission(name, educt1, product1, product2, rate, distance, weight1, weight2);
                } else {
                    // enzymatic
                    std::tie(product1, product2, productDistance) = treatSideWithPlus(rhs, plusPos, false);
                    if (educt2 != product2) {
                        throw std::invalid_argument(fmt::format(
                                R"(In enzymatic reaction, educt2 ("{}") and product2 ("{}") have to be equal)",
                                educt2, product2));
                    }
                    auto distance = static_cast<scalar>(std::stod(eductDistance));
                    return addEnzymatic(name, educt2, educt1, product1, rate, distance);
                }
            }
        }
    }
}

namespace {
struct FindId {
    explicit FindId(std::string name, ReactionId &id) : name(std::move(name)), id(id) {}

    std::string name {""};
    std::reference_wrapper<ReactionId> id;
    std::shared_ptr<bool> found = std::make_shared<bool>(false);

    template<typename T>
    void operator()(const T &reaction) {
        if (reaction->name() == name) {
            *found = true;
            id.get() = reaction->id();
        }
    }
};

struct FindName {
    explicit FindName(ReactionId id, std::string &name) : id(id), name(name) {}

    ReactionId id {0};
    std::reference_wrapper<std::string> name;
    std::shared_ptr<bool> found = std::make_shared<bool>(false);

    template<typename T>
    void operator()(const T &reaction) {
        if (reaction->id() == id) {
            *found = true;
            name.get() = reaction->name();
        };
    }
};

struct FindReactionById {
    explicit FindReactionById(ReactionId id) : id(id) {}

    bool found() const {
        return _found;
    }

    template<typename T>
    void operator()(const T &r) {
        if (r->id() == id) {
            if (!_found) {
                _found = true;
                assign(r);
            } else {
                throw std::runtime_error(
                        fmt::format("reaction was already found, there cannot exist two with the same id ({})", id));
            }
        };
    }

    const Reaction* result() {
        if (_found) {
            return reaction;
        } else {
            throw std::runtime_error(fmt::format("reaction with id {} was not found", id));
        }
    }

private:
    void assign(const Reaction* r2) {
        reaction = r2;
    }

    void assign(const std::shared_ptr<const Reaction> &r2) {
        reaction = r2.get();
    }

    const ReactionId id;
    const Reaction* reaction = nullptr;
    bool _found = false;
};
}

bool ReactionRegistry::reactionNameExists(const std::string &name) const {
    ReactionId id;
    FindId findId(name, id);
    readdy::util::collections::for_each_value_ref(_ownO1Reactions, findId);
    readdy::util::collections::for_each_value_ref(_ownO2Reactions, findId);
    return *(findId.found);
}

std::string ReactionRegistry::nameOf(ReactionId id) const {
    std::string name;
    FindName findName(id, name);
    readdy::util::collections::for_each_value_ref(_ownO1Reactions, findName);
    readdy::util::collections::for_each_value_ref(_ownO2Reactions, findName);
    if (*(findName.found)) {
        return name;
    } else {
        throw std::runtime_error(fmt::format("no reaction with id {} exists", id));
    }
}

ReactionId ReactionRegistry::idOf(const std::string &name) const {
    ReactionId id;
    FindId findId(name, id);
    readdy::util::collections::for_each_value_ref(_ownO1Reactions, findId);
    readdy::util::collections::for_each_value_ref(_ownO2Reactions, findId);
    if (*(findId.found)) {
        return findId.id;
    } else {
        throw std::runtime_error(fmt::format("no reaction with name {} exists", name));
    }
}

const Reaction* ReactionRegistry::byId(ReactionId id) const {
    FindReactionById findReaction(id);
    readdy::util::collections::for_each_value_ref(_ownO1Reactions, findReaction);
    readdy::util::collections::for_each_value_ref(_ownO2Reactions, findReaction);
    if (findReaction.found()) {
        return findReaction.result();
    } else {
        throw std::runtime_error(fmt::format("no reaction with id {} exists", id));
    }
}

}