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
 * @file SCPUReactionUtils.h
 * @brief << brief description >>
 * @author clonker
 * @author chrisfroe
 * @date 31.01.17
 * @copyright BSD-3
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
            bcs::fixPosition(p.pos(), box, pbc);
            newEntries.emplace_back(p);

            entry1.type = reaction->products()[0];
            entry1.pos += reaction->weight1() * distance * n3;
            entry1.id = readdy::model::Particle::nextId();
            bcs::fixPosition(entry1.pos, box, pbc);
            if(record) {
                record->products[0] = entry1.id;
                record->products[1] = p.id();
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
