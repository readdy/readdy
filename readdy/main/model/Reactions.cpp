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
 * @file Reactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.10.16
 */

#include "readdy/model/reactions/Reaction.h"
#include "readdy/model/reactions/ReactionRecord.h"

namespace readdy {
namespace model {
namespace reactions {

template<> short Reaction<1>::counter = 0;
template<> short Reaction<2>::counter = 0;

std::ostream& operator<<(std::ostream& os, const ReactionType& reactionType) {
    switch (reactionType) {
        case ReactionType::Decay: os << "Decay"; break;
        case ReactionType::Conversion: os << "Conversion"; break;
        case ReactionType::Fusion: os << "Fusion"; break;
        case ReactionType::Fission: os << "Fission"; break;
        case ReactionType::Enzymatic: os << "Enzymatic"; break;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const ReactionRecord& record) {
    ReactionType type = ReactionType(record.type);
    os << "ReactionRecord[type: " << type;
    switch (type) {
        case ReactionType::Decay:{
            os << " educt: " << record.educts[0];
            break;
        }
        case ReactionType::Conversion: {
            os << " educt: " << record.educts[0] << " product: " << record.products[0];
            break;
        }
        case ReactionType::Fusion: {
            os << " educts: " << record.educts[0] << "," << record.educts[1] << " product: " << record.products[0];
            break;
        }
        case ReactionType::Fission: {
            os << " educt: " << record.educts[0] << " products: " << record.products[0] << "," << record.products[1];
            break;
        }
        case ReactionType::Enzymatic: {
            os << " educts: " << record.educts[0] << "," << record.educts[1];
            os << " products: " << record.products[0] << "," << record.products[1];
            break;
        }
    }
    os << " location: " << record.where;
    return os;
}

}
}
}