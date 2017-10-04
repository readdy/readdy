/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file ReadableReactionRecord.h
 * @brief << brief description >>
 * @author clonker
 * @date 04.10.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <sstream>
#include <readdy/model/Particle.h>

namespace rpy {

namespace detail {

inline std::string nameOf(const readdy::model::ioutils::ReactionInfo &info) {
    return info.name;
}

inline std::string nameOf(readdy::model::reactions::Reaction<1> *reaction) {
    return reaction->name();
}

inline std::string nameOf(const readdy::model::reactions::Reaction<2> *reaction) {
    return reaction->name();
}

}

struct ReadableReactionRecord {
    ReadableReactionRecord() : type(""), reaction_label(""), educts(), products(), where() {}

    std::string type;
    std::string reaction_label;
    std::vector<readdy::model::Particle::id_type> educts;
    std::vector<readdy::model::Particle::id_type> products;
    readdy::Vec3::data_arr where;
};

template<typename ReactionsOrder1, typename ReactionsOrder2>
inline ReadableReactionRecord convert(const readdy::model::reactions::ReactionRecord &reaction,
                               const ReactionsOrder1& reactionInfoOrder1, const ReactionsOrder2 &reactionInfoOrder2) {
    rpy::ReadableReactionRecord rrr {};
    auto tt = readdy::model::reactions::ReactionType(reaction.type);
    switch(tt) {
        case readdy::model::reactions::ReactionType::Conversion:{
            rrr.educts = {reaction.educts[0]};
            rrr.products = {reaction.products[0]};
            rrr.where = reaction.where.data;
            rrr.type = "conversion";
            rrr.reaction_label = detail::nameOf(reactionInfoOrder1.at(reaction.reactionIndex));
            break;
        };
        case readdy::model::reactions::ReactionType::Fusion: {
            rrr.educts = {reaction.educts[0], reaction.educts[1]};
            rrr.products = {reaction.products[0]};
            rrr.where = reaction.where.data;
            rrr.type = "fusion";
            rrr.reaction_label = detail::nameOf(reactionInfoOrder2.at(reaction.reactionIndex));
            break;
        };
        case readdy::model::reactions::ReactionType::Fission: {
            rrr.educts = {reaction.educts[0]};
            rrr.products = {reaction.products[0], reaction.products[1]};
            rrr.where = reaction.where.data;
            rrr.type = "fission";
            rrr.reaction_label = detail::nameOf(reactionInfoOrder1.at(reaction.reactionIndex));
            break;
        };
        case readdy::model::reactions::ReactionType::Enzymatic: {
            rrr.educts = {reaction.educts[0], reaction.educts[1]};
            rrr.products = {reaction.products[0], reaction.products[1]};
            rrr.where = reaction.where.data;
            rrr.type = "enzymatic";
            rrr.reaction_label = detail::nameOf(reactionInfoOrder2.at(reaction.reactionIndex));
            break;
        };
        case readdy::model::reactions::ReactionType::Decay: {
            rrr.educts = {reaction.educts[0]};
            rrr.products = {};
            rrr.where = reaction.where.data;
            rrr.type = "decay";
            rrr.reaction_label = detail::nameOf(reactionInfoOrder1.at(reaction.reactionIndex));
            break;
        };
    }

    return rrr;
}

inline std::string repr(const ReadableReactionRecord &rrr) {
    std::stringstream result;

    result << "Reaction[type=" << rrr.type << ", label="<< rrr.reaction_label << ", educts=[";
    if(rrr.educts.size() == 1) {
        result << rrr.educts.at(0) << "]";
    } else {
        result << rrr.educts.at(0) << ", " << rrr.educts.at(1) << "]";
    }
    result << ", products=[";
    if(rrr.products.size() == 1) {
        result << rrr.products.at(0) << "]";
    } else if (rrr.products.size() == 2){
        result << rrr.products.at(0) << ", " << rrr.products.at(1) << "]";
    }
    result << ", position=(" << rrr.where.at(0) << ", " << rrr.where.at(1) << ", " << rrr.where.at(2) << ")]";
    return result.str();
}

}