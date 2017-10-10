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

struct ReadableReactionRecord {
    ReadableReactionRecord();

    ReadableReactionRecord(const ReadableReactionRecord &) = default;

    ReadableReactionRecord &operator=(const ReadableReactionRecord &) = default;

    ReadableReactionRecord(ReadableReactionRecord &&) = default;

    ReadableReactionRecord &operator=(ReadableReactionRecord &&) = default;

    ~ReadableReactionRecord() = default;

    std::string type;
    std::string reaction_label;
    std::vector<readdy::model::Particle::id_type> educts;
    std::vector<readdy::model::Particle::id_type> products;
    readdy::Vec3::data_arr where;

    friend std::ostream &operator<<(std::ostream &os, const ReadableReactionRecord &rrr);
};

ReadableReactionRecord convert(const readdy::model::reactions::ReactionRecord &reaction, const std::string &name);

std::string repr(const ReadableReactionRecord &rrr);

}
