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
 * @file SCPUReactionUtils.h
 * @brief << brief description >>
 * @author clonker
 * @author chrisfroe
 * @date 31.01.17
 * @copyright GNU Lesser General Public License v3.0
 */
#pragma once

#include <readdy/common/boundary_condition_operations.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
namespace reactions {

using scpu_model = readdy::kernel::scpu::SCPUStateModel;
using scpu_data = readdy::kernel::scpu::model::SCPUParticleData;
using reaction_type = readdy::model::reactions::ReactionType;
using context = readdy::model::Context;
using fix_pos = context::fix_pos_fun;
using reaction_record = readdy::model::reactions::ReactionRecord;

template<typename Reaction>
void performReaction(
        scpu_data &data, scpu_data::entry_index idx1, scpu_data::entry_index idx2, scpu_data::new_entries &newEntries,
        std::vector<scpu_data::entry_index> &decayedEntries, Reaction *reaction, const readdy::model::Context& context, reaction_record* record) {
    //const auto& fixPos = context.fixPositionFun();
    //const auto& shortestDifferenceFun = context.shortestDifferenceFun();
    const auto &box = context.boxSize().data();
    const auto &pbc = context.periodicBoundaryConditions().data();
    auto& entry1 = data.entry_at(idx1);
    auto& entry2 = data.entry_at(idx2);
    if(record) {
        record->type = static_cast<int>(reaction->type());
        record->where = (entry1.position() + entry2.position()) / 2.;
        bcs::fixPosition(record->where, box, pbc);
        record->educts[0] = entry1.id;
        record->educts[1] = entry2.id;
        record->types_from[0] = entry1.type;
        record->types_from[1] = entry2.type;
    }
    switch (reaction->type()) {
        case reaction_type::Decay: {
            decayedEntries.push_back(idx1);
            break;
        }
        case reaction_type::Conversion: {
            entry1.type = reaction->products()[0];
            entry1.id = readdy::model::Particle::nextId();
            if(record) record->products[0] = entry1.id;
            break;
        }
        case reaction_type::Enzymatic: {
            if (entry1.type == reaction->educts()[1]) {
                // p1 is the catalyst
                entry2.type = reaction->products()[0];
                entry2.id = readdy::model::Particle::nextId();
            } else {
                // p2 is the catalyst
                entry1.type = reaction->products()[0];
                entry1.id = readdy::model::Particle::nextId();
            }
            if(record) {
                record->products[0] = entry1.id;
                record->products[1] = entry2.id;
            }
            break;
        }
        case reaction_type::Fission: {
            auto n3 = readdy::model::rnd::normal3<readdy::scalar>(0, 1);
            n3 /= std::sqrt(n3 * n3);

            const auto distance =
                    reaction->productDistance() * std::cbrt(readdy::model::rnd::uniform_real<readdy::scalar>(0, 1));
            readdy::model::Particle p(entry1.position() - reaction->weight2() * distance * n3,
                                      reaction->products()[1]);
            bcs::fixPosition(p.getPos(), box, pbc);
            newEntries.emplace_back(p);

            entry1.type = reaction->products()[0];
            entry1.pos += reaction->weight1() * distance * n3;
            entry1.id = readdy::model::Particle::nextId();
            bcs::fixPosition(entry1.pos, box, pbc);
            if(record) {
                record->products[0] = entry1.id;
                record->products[1] = p.getId();
            }
            break;
        }
        case reaction_type::Fusion: {
            const auto e1Pos = data.entry_at(idx1).pos;
            const auto e2Pos = data.entry_at(idx2).pos;
            const auto difference = bcs::shortestDifference(e1Pos, e2Pos, box, pbc);
            if (reaction->educts()[0] == entry1.type) {
                entry1.pos += reaction->weight1() * difference;
            } else {
                entry1.pos += reaction->weight2() * difference;
            }
            bcs::fixPosition(entry1.pos, box, pbc);
            entry1.type = reaction->products()[0];
            entry1.id = readdy::model::Particle::nextId();
            decayedEntries.push_back(idx2);
            if(record) record->products[0] = entry1.id;
            break;
        }
    }
}

}
}
}
}
}
