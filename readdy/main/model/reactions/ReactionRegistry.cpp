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
 * @author chrisfroe
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
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

namespace readdy {
namespace model {
namespace reactions {

ReactionRegistry::reaction_id ReactionRegistry::emplaceReaction(const std::shared_ptr<Reaction> &reaction) {
    if (reactionNameExists(reaction->name())) {
        throw std::invalid_argument(fmt::format("A reaction with the name {} exists already", reaction->name()));
    }
    reaction_id id = reaction->id();
    if (reaction->nEducts() == 1) {
        auto t = reaction->educts()[0];
        if (one_educt_registry_internal.find(t) == one_educt_registry_internal.end()) {
            one_educt_registry_internal.emplace(t, rea_ptr_vec());
        }
        one_educt_registry_internal[t].push_back(reaction);
        _n_order1 += 1;
    } else {
        const auto pp = std::make_tuple(reaction->educts()[0], reaction->educts()[1]);
        if (two_educts_registry_internal.find(pp) == two_educts_registry_internal.end()) {
            two_educts_registry_internal.emplace(pp, rea_ptr_vec());
        }
        two_educts_registry_internal[pp].push_back(reaction);
        _n_order2 += 1;
    }
    return id;
}

void ReactionRegistry::configure() {
    namespace coll = readdy::util::collections;
    using pair = readdy::util::particle_type_pair;
    using reaction1ptr = std::shared_ptr<reaction>;
    using reaction2ptr = std::shared_ptr<reaction>;

    one_educt_registry.clear();
    two_educts_registry.clear();
    _reaction_o2_types.clear();

    coll::for_each_value(one_educt_registry_internal,
                         [&](const particle_type_type type, const reaction1ptr &ptr) {
                             (one_educt_registry)[type].push_back(ptr.get());
                         });
    coll::for_each_value(two_educts_registry_internal, [&](const pair &type, const reaction2ptr &r) {
        (two_educts_registry)[type].push_back(r.get());
        _reaction_o2_types.emplace(std::get<0>(type));
        _reaction_o2_types.emplace(std::get<1>(type));
    });
    coll::for_each_value(one_educt_registry_external,
                         [&](const particle_type_type type, reaction *ptr) {
                             (one_educt_registry)[type].push_back(ptr);
                         });
    coll::for_each_value(two_educts_registry_external, [&](const pair &type, reaction *r) {
        (two_educts_registry)[type].push_back(r);
        _reaction_o2_types.emplace(std::get<0>(type));
        _reaction_o2_types.emplace(std::get<1>(type));
    });

}

std::string ReactionRegistry::describe() const {
    namespace rus = readdy::util::str;
    std::string description;
    auto nameOf = [&](particle_type_type t) {return _types.get().nameOf(t);};

    if (!one_educt_registry.empty()) {
        description += fmt::format(" - unimolecular reactions:{}", rus::newline);
        for (const auto &entry : one_educt_registry) {
            for (const auto reaction : entry.second) {
                const auto &educts = reaction->educts();
                const auto &products = reaction->products();

                switch (reaction->type()) {
                    case ReactionType::Conversion: {
                        description += fmt::format("     * Conversion {} -> {} with a rate of {}{}",
                                                   nameOf(educts[0]), nameOf(products[0]), reaction->rate(),
                                                   rus::newline);
                        break;
                    }
                    case ReactionType::Fission: {
                        description += fmt::format("     * Fission {} -> {} + {} with a rate of {}, a product distance of {}, and weights {} and {}{}",
                                                   nameOf(educts[0]), nameOf(products[0]), nameOf(products[1]),
                                                   reaction->rate(), reaction->productDistance(),
                                                   reaction->weight1(), reaction->weight2(), rus::newline);
                        break;
                    }
                    case ReactionType::Decay: {
                        description += fmt::format("     * Decay {} -> ø with a rate of {}{}",
                                                   nameOf(educts[0]), reaction->rate(), rus::newline);
                        break;
                    }
                    default: {
                        throw std::logic_error(fmt::format("unimolecular reaction in registry that was of type {}", 
                                                           reaction->type()));
                    }
                }
            }
        }
    }
    if (!two_educts_registry.empty()) {
        description += fmt::format(" - bimolecular reactions:{}", rus::newline);
        for (const auto &entry : two_educts_registry) {
            for (const auto reaction : entry.second) {
                const auto &educts = reaction->educts();
                const auto &products = reaction->products();

                switch(reaction->type()) {
                    case ReactionType::Enzymatic: {
                        auto enzymatic = dynamic_cast<const Enzymatic*>(reaction);
                        description += fmt::format("     * Enzymatic {} + {} -> {} + {} with a rate of {} and an educt distance of {}{}",
                                                   nameOf(enzymatic->getFrom()), nameOf(enzymatic->getCatalyst()),
                                                   nameOf(enzymatic->getTo()), nameOf(enzymatic->getCatalyst()),
                                                   enzymatic->rate(), enzymatic->eductDistance(), rus::newline);
                        break;
                    }
                    case ReactionType::Fusion: {
                        auto fusion = dynamic_cast<const Fusion*>(reaction);
                        description += fmt::format("     * Fusion {} + {} -> {} with a rate of {}, an educt distance of {}, and weights {} and {}{}",
                                                   nameOf(fusion->getFrom1()), nameOf(fusion->getFrom2()),
                                                   nameOf(fusion->getTo()), fusion->rate(), fusion->eductDistance(),
                                                   fusion->weight1(), fusion->weight2(), rus::newline);
                        break;
                    }
                    default: {
                        throw std::logic_error(fmt::format("bimolecular reaction in registry that was of type {}",
                                                           reaction->type()));
                    }
                }
            }
        }
    }
    return description;
}

const short ReactionRegistry::addExternal(reaction *r) {
    if (r->nEducts() == 1) {
        one_educt_registry_external[r->educts()[0]].push_back(r);
        _n_order1 += 1;
        return r->id();
    }
    two_educts_registry_external[std::tie(r->educts()[0], r->educts()[1])].push_back(r);
    _n_order2 += 1;
    return r->id();
}

ReactionRegistry::reaction_id ReactionRegistry::add(const std::string &descriptor, scalar rate) {
    namespace mutil = readdy::model::util;
    namespace rutil = readdy::util;
    log::trace("begin parsing \"{}\"", descriptor);
    auto arrowPos = descriptor.find(mutil::arrow);
    if (arrowPos == descriptor.npos) {
        throw std::invalid_argument(fmt::format(
                "the descriptor must contain an arrow (\"{}\") to indicate lhs and rhs.", mutil::arrow
        ));
    }
    if (descriptor.find(mutil::arrow, arrowPos + 1) != descriptor.npos) {
        throw std::invalid_argument(fmt::format(
                "the descriptor must not contain more than one arrow (\"{}\").", mutil::arrow
        ));
    }
    auto lhs = descriptor.substr(0, arrowPos);
    auto rhs = descriptor.substr(arrowPos + std::strlen(mutil::arrow), descriptor.npos);

    rutil::str::trim(lhs);
    rutil::str::trim(rhs);

    std::string name;
    {
        auto colonPos = lhs.find(':');
        if (colonPos == lhs.npos) {
            throw std::invalid_argument("The descriptor did not contain a colon ':' to specify the end of the name.");
        }
        name = rutil::str::trim_copy(lhs.substr(0, colonPos));
        lhs = rutil::str::trim_copy(lhs.substr(colonPos + 1, lhs.npos));
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
            if (commaPos == weightsString.npos) {
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
            if (closingParentheses == s.npos) {
                throw std::invalid_argument(fmt::format("The side (\"{}\") did not contain a ')'.", s));
            }
            pType2 = rutil::str::trim_copy(s.substr(closingParentheses + 1, s.npos));
            reactionRadius = getReactionRadius(s);
        } else {
            pType2 = rutil::str::trim_copy(s.substr(plusPos + 1, s.npos));
        }
        return std::make_tuple(pType1, pType2, reactionRadius);
    };

    std::string w1, w2;
    scalar weight1{0.5}, weight2{0.5};
    {
        auto bracketPos = rhs.find('[');
        if (bracketPos != rhs.npos) {
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
        if (plusPos == lhs.npos) {
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
            if (plusPos == rhs.npos) {
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
    explicit FindId(std::string name, Reaction::reaction_id &id) : name(std::move(name)), id(id) {}

    std::string name {""};
    std::reference_wrapper<Reaction::reaction_id> id;
    std::shared_ptr<bool> found = std::make_shared<bool>(false);

    template<typename T>
    void operator()(const T &reaction) {
        if (reaction->name() == name) {
            *found = true;
            id.get() = reaction->id();
        };
    }
};

struct FindName {
    explicit FindName(Reaction::reaction_id id, std::string &name) : id(id), name(name) {}

    Reaction::reaction_id id {0};
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
}

bool ReactionRegistry::reactionNameExists(const std::string &name) const {
    Reaction::reaction_id id;
    FindId findId(name, id);
    readdy::util::collections::for_each_value_ref(one_educt_registry_internal, findId);
    readdy::util::collections::for_each_value_ref(one_educt_registry_external, findId);
    readdy::util::collections::for_each_value_ref(two_educts_registry_internal, findId);
    readdy::util::collections::for_each_value_ref(two_educts_registry_external, findId);
    return *(findId.found);
}

std::string ReactionRegistry::nameOf(ReactionRegistry::reaction_id id) const {
    std::string name;
    FindName findName(id, name);
    readdy::util::collections::for_each_value_ref(one_educt_registry_internal, findName);
    readdy::util::collections::for_each_value_ref(one_educt_registry_external, findName);
    readdy::util::collections::for_each_value_ref(two_educts_registry_internal, findName);
    readdy::util::collections::for_each_value_ref(two_educts_registry_external, findName);
    if (*(findName.found)) {
        return name;
    } else {
        throw std::runtime_error(fmt::format("no reaction with id {} exists", id));
    }
}

ReactionRegistry::reaction_id ReactionRegistry::idOf(const std::string &name) const {
    Reaction::reaction_id id;
    FindId findId(name, id);
    readdy::util::collections::for_each_value_ref(one_educt_registry_internal, findId);
    readdy::util::collections::for_each_value_ref(one_educt_registry_external, findId);
    readdy::util::collections::for_each_value_ref(two_educts_registry_internal, findId);
    readdy::util::collections::for_each_value_ref(two_educts_registry_external, findId);
    if (*(findId.found)) {
        return findId.id;
    } else {
        throw std::runtime_error(fmt::format("no reaction with name {} exists", name));
    }
}

}
}
}