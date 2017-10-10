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
 * @file TopologyFusionReaction.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/reactions/SpatialTopologyReaction.h>

#include <regex>
#include <readdy/model/topologies/TopologyRegistry.h>
#include <readdy/model/Utils.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {

SpatialTopologyReaction::SpatialTopologyReaction(std::string name, readdy::util::particle_type_pair types,
                                                 topology_type_pair top_types,
                                                 readdy::util::particle_type_pair types_to, topology_type_pair top_types_to,
                                                 scalar rate, scalar radius, STRMode mode)
        : _name(std::move(name)), _types(std::move(types)), _types_to(std::move(types_to)), _rate(rate),
          _radius(radius), _mode(mode), _top_types(std::move(top_types)), _top_types_to(std::move(top_types_to)) {}

const std::string &SpatialTopologyReaction::name() const {
    return _name;
}

const particle_type_type SpatialTopologyReaction::type1() const {
    return std::get<0>(_types);
}

const particle_type_type SpatialTopologyReaction::type2() const {
    return std::get<1>(_types);
}

const scalar SpatialTopologyReaction::rate() const {
    return _rate;
}

const scalar SpatialTopologyReaction::radius() const {
    return _radius;
}

const readdy::util::particle_type_pair &SpatialTopologyReaction::types() const {
    return _types;
}

const particle_type_type SpatialTopologyReaction::type_to1() const {
    return std::get<0>(_types_to);
}

const particle_type_type SpatialTopologyReaction::type_to2() const {
    return std::get<1>(_types_to);
}

const readdy::util::particle_type_pair &SpatialTopologyReaction::types_to() const {
    return _types_to;
}

const bool SpatialTopologyReaction::allow_self_connection() const {
    return _mode == STRMode::TT_FUSION_ALLOW_SELF;
}

const topology_type_type SpatialTopologyReaction::top_type1() const {
    return std::get<0>(_top_types);
}

const topology_type_type SpatialTopologyReaction::top_type2() const {
    return std::get<1>(_top_types);
}

const topology_type_type SpatialTopologyReaction::top_type_to1() const {
    return std::get<0>(_top_types_to);
}

const topology_type_type SpatialTopologyReaction::top_type_to2() const {
    return std::get<1>(_top_types_to);
}

bool SpatialTopologyReaction::is_topology_particle_reaction() const {
    return top_type2() == topology_type_empty;
}

bool SpatialTopologyReaction::is_topology_topology_reaction() const {
    return !is_topology_particle_reaction();
}

const bool SpatialTopologyReaction::is_enzymatic() const {
    return _mode == STRMode::TT_ENZYMATIC || _mode == STRMode::TP_ENZYMATIC;
}

const bool SpatialTopologyReaction::is_fusion() const {
    return _mode == STRMode::TT_FUSION || _mode == STRMode::TT_FUSION_ALLOW_SELF || _mode == STRMode::TP_FUSION;
}

const STRMode &SpatialTopologyReaction::mode() const {
    return _mode;
}

SpatialTopologyReaction STRParser::parse(const std::string &descriptor, scalar rate, scalar radius) const {
    namespace mutil = readdy::model::util;
    namespace rutil = readdy::util;
    SpatialTopologyReaction reaction;
    reaction._rate = rate;
    reaction._radius = radius;

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
    reaction._name = name;

    static std::regex particleTypeRegex(R"(\(([^\)]*)\))");
    static std::regex topologyTypeRegex(R"([^(]*)");

    static auto getTop = [](const std::string &s) {
        std::smatch topMatch;
        if (std::regex_search(s, topMatch, topologyTypeRegex)) {
            return rutil::str::trim_copy(topMatch.str());
        }
        throw std::invalid_argument(fmt::format("The term \"{}\" did not contain a topology type.", s));
    };
    static auto getParticleType = [](const std::string &s) {
        std::smatch ptMatch;
        if (std::regex_search(s, ptMatch, particleTypeRegex)) {
            auto pt = rutil::str::trim_copy(ptMatch.str());
            return rutil::str::trim_copy(pt.substr(1, pt.size() - 2));
        }
        throw std::invalid_argument(fmt::format("The term \"{}\" did not contain a particle type.", s));
    };

    static auto treatTerm = [](const std::string &s) {
        return std::make_tuple(getParticleType(s), getTop(s));
    };

    static auto treatSide = [](const std::string &s) {
        auto plusPos = s.find('+');
        if (plusPos == s.npos) {
            throw std::invalid_argument("The left hand side of the topology reaction did not contain a '+'.");
        }
        auto educt1 = rutil::str::trim_copy(s.substr(0, plusPos));
        auto educt2 = rutil::str::trim_copy(s.substr(plusPos + 1, s.npos));

        std::string t1, t2, p1, p2;
        std::tie(p1, t1) = treatTerm(educt1);
        std::tie(p2, t2) = treatTerm(educt2);

        return std::make_tuple(p1, t1, p2, t2);
    };

    std::string lhs_p1, lhs_t1, lhs_p2, lhs_t2;
    {
        std::tie(lhs_p1, lhs_t1, lhs_p2, lhs_t2) = treatSide(lhs);
    }

    std::string rhs_p1, rhs_t1, rhs_p2, rhs_t2;
    bool rhs_fusion {false};
    {
        auto plusPos = rhs.find('+');
        rhs_fusion = plusPos == rhs.npos;
        if (plusPos == rhs.npos) {
            // fusion type
            std::string fuse;
            std::tie(fuse, rhs_t1) = treatTerm(rhs);
            rhs_t2 = "";

            auto separatorPos = fuse.find(mutil::bond);
            if (separatorPos == fuse.npos) {
                throw std::invalid_argument(fmt::format(
                        "The right-hand side was of fusion type but there was no bond \"{}\" defined.", mutil::bond
                ));
            }

            rhs_p1 = rutil::str::trim_copy(fuse.substr(0, separatorPos));
            rhs_p2 = rutil::str::trim_copy(fuse.substr(separatorPos + std::strlen(mutil::bond), fuse.npos));
        } else {
            // enzymatic type
            std::tie(rhs_p1, rhs_t1, rhs_p2, rhs_t2) = treatSide(rhs);
        }

        log::trace(R"(got lhs with toplogies "{}" and "{}", particle types "{}" and "{}")", lhs_t1, lhs_t2, lhs_p1,
                   lhs_p2);
        log::trace(R"(got rhs with topologies "{}" and "{}", particles "{}" and "{}")", rhs_t1, rhs_t2, rhs_p1, rhs_p2);
    }

    const auto &particle_types = _topology_registry.get().particleTypeRegistry();

    if(lhs_t2.empty()) {
        // we are in the topology-particle case
        reaction._top_types = std::make_tuple(_topology_registry.get().idOf(lhs_t1), topology_type_empty);
        reaction._types = std::make_tuple(particle_types.idOf(lhs_p1), particle_types.idOf(lhs_p2));
        reaction._types_to = std::make_tuple(particle_types.idOf(rhs_p1), particle_types.idOf(rhs_p2));
        reaction._top_types_to = std::make_tuple(_topology_registry.get().idOf(rhs_t1), topology_type_empty);
        if(rhs_fusion) {
            // we are in the fusion case
            reaction._mode = STRMode::TP_FUSION;
        } else {
            // we are in the enzymatic case
            reaction._mode = STRMode::TP_ENZYMATIC;
        }
    } else {
        // we are in the topology-topology case
        reaction._top_types = std::make_tuple(_topology_registry.get().idOf(lhs_t1),
                                              _topology_registry.get().idOf(lhs_t2));
        reaction._types = std::make_tuple(particle_types.idOf(lhs_p1), particle_types.idOf(lhs_p2));
        reaction._types_to = std::make_tuple(particle_types.idOf(rhs_p1), particle_types.idOf(rhs_p2));
        if(rhs_fusion) {
            // we are in the fusion case
            if(rhs.find("[self=true]") != rhs.npos) {
                reaction._mode = STRMode::TT_FUSION_ALLOW_SELF; // allow self?
            } else {
                reaction._mode = STRMode::TT_FUSION;
            }
            reaction._top_types_to = std::make_tuple(_topology_registry.get().idOf(rhs_t1), topology_type_empty);
        } else {
            // we are in the enzymatic case
            reaction._mode = STRMode::TT_ENZYMATIC;
            reaction._top_types_to = std::make_tuple(_topology_registry.get().idOf(rhs_t1),
                                                     _topology_registry.get().idOf(rhs_t2));
        }
    }

    return reaction;
}

STRParser::STRParser(const TopologyRegistry &registry) : _topology_registry(registry) {}

}
}
}
}