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
 * @file ReadableReactionRecord.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.10.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/reactions/ReactionRecord.h>
#include "ReadableReactionRecord.h"

namespace rpy {

ReadableReactionRecord convert(const readdy::model::reactions::ReactionRecord &reaction, const std::string &name) {
    rpy::ReadableReactionRecord rrr{};
    rrr.where = reaction.where.data;
    rrr.reaction_label = name;
    auto tt = readdy::model::reactions::ReactionType(reaction.type);
    switch (tt) {
        case readdy::model::reactions::ReactionType::Conversion: {
            rrr.educts = {reaction.educts[0]};
            rrr.products = {reaction.products[0]};
            rrr.type = "conversion";
            break;
        };
        case readdy::model::reactions::ReactionType::Fusion: {
            rrr.educts = {reaction.educts[0], reaction.educts[1]};
            rrr.products = {reaction.products[0]};
            rrr.type = "fusion";
            break;
        };
        case readdy::model::reactions::ReactionType::Fission: {
            rrr.educts = {reaction.educts[0]};
            rrr.products = {reaction.products[0], reaction.products[1]};
            rrr.type = "fission";
            break;
        };
        case readdy::model::reactions::ReactionType::Enzymatic: {
            rrr.educts = {reaction.educts[0], reaction.educts[1]};
            rrr.products = {reaction.products[0], reaction.products[1]};
            rrr.type = "enzymatic";
            break;
        };
        case readdy::model::reactions::ReactionType::Decay: {
            rrr.educts = {reaction.educts[0]};
            rrr.products = {};
            rrr.type = "decay";
            break;
        };
    }
    return rrr;
}

std::string repr(const ReadableReactionRecord &rrr) {
    std::stringstream result;

    result << "Reaction[type=" << rrr.type << ", label=" << rrr.reaction_label << ", educts=[";
    if (rrr.educts.size() == 1) {
        result << rrr.educts.at(0) << "]";
    } else {
        result << rrr.educts.at(0) << ", " << rrr.educts.at(1) << "]";
    }
    result << ", products=[";
    if (rrr.products.size() == 1) {
        result << rrr.products.at(0) << "]";
    } else if (rrr.products.size() == 2) {
        result << rrr.products.at(0) << ", " << rrr.products.at(1) << "]";
    }
    result << ", position=(" << rrr.where.at(0) << ", " << rrr.where.at(1) << ", " << rrr.where.at(2) << ")]";
    return result.str();
}

std::ostream &operator<<(std::ostream &os, const ReadableReactionRecord &rrr) {
    os << repr(rrr);
    return os;
}

ReadableReactionRecord::ReadableReactionRecord()
        : type(""), reaction_label(""), educts(), products(), where({0, 0, 0}) {}

}
