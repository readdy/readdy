/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file ObservableData.h
 * @brief << brief description >>
 * @author clonker
 * @date 1/17/18
 */


#pragma once

#include <readdy/common/common.h>
#include <readdy/model/reactions/ReactionRecord.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {

struct ObservableData {

    using reaction_counts_map = readdy::model::reactions::reaction_counts_map;

    scalar energy = 0;
    std::vector<readdy::model::reactions::ReactionRecord> reactionRecords{};
    reaction_counts_map reactionCounts {};
    Matrix33 virial {};
};

}
}
}
}
