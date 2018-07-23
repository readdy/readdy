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
 * @file ReadableReactionRecord.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.10.17
 * @copyright GPL-3
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
