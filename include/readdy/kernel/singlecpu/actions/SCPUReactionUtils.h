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
 * @date 31.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#ifndef READDY_MAIN_SCPUREACTIONUTILS_H
#define READDY_MAIN_SCPUREACTIONUTILS_H
namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
namespace reactions {

using model_t = readdy::kernel::scpu::SCPUStateModel;
using data_t = readdy::kernel::scpu::model::SCPUParticleData;
using reaction_type = readdy::model::reactions::ReactionType;
using ctx_t = readdy::model::KernelContext;
using fix_pos = ctx_t::fix_pos_fun;
using record_t = readdy::model::reactions::ReactionRecord;

template<typename Reaction>
void performReaction(
        data_t &data, data_t::index_t idx1, data_t::index_t idx2, data_t::entries_update_t &newEntries,
        std::vector<data_t::index_t> &decayedEntries, Reaction *reaction, const fix_pos &fixPos, record_t* record) {
    auto& entry1 = data.entry_at(idx1);
    auto& entry2 = data.entry_at(idx2);
    if(record) {
        record->type = static_cast<int>(reaction->getType());
        record->where = (entry1.position() + entry2.position()) / 2.;
        fixPos(record->where);
        record->educts[0] = entry1.id;
        record->educts[1] = entry2.id;
        record->types_from[0] = entry1.type;
        record->types_from[1] = entry2.type;
    }
    switch (reaction->getType()) {
        case reaction_type::Decay: {
            decayedEntries.push_back(idx1);
            break;
        }
        case reaction_type::Conversion: {
            entry1.type = reaction->getProducts()[0];
            entry1.id = readdy::model::Particle::nextId();
            if(record) record->products[0] = entry1.id;
            break;
        }
        case reaction_type::Enzymatic: {
            if (entry1.type == reaction->getEducts()[1]) {
                // p1 is the catalyst
                entry2.type = reaction->getProducts()[0];
                entry2.id = readdy::model::Particle::nextId();
            } else {
                // p2 is the catalyst
                entry1.type = reaction->getProducts()[0];
                entry1.id = readdy::model::Particle::nextId();
            }
            if(record) {
                record->products[0] = entry1.id;
                record->products[1] = entry2.id;
            }
            break;
        }
        case reaction_type::Fission: {
            auto n3 = readdy::model::rnd::normal3(0, 1);
            n3 /= std::sqrt(n3 * n3);

            readdy::model::Particle p(entry1.position() - reaction->getWeight2() * reaction->getProductDistance() * n3,
                                      reaction->getProducts()[1]);
            newEntries.push_back({p});

            entry1.type = reaction->getProducts()[0];
            entry1.pos += reaction->getWeight1() * reaction->getProductDistance() * n3;
            entry1.id = readdy::model::Particle::nextId();
            fixPos(entry1.pos);
            if(record) {
                record->products[0] = entry1.id;
                record->products[1] = p.getId();
            }
            break;
        }
        case reaction_type::Fusion: {
            const auto e1Pos = data.entry_at(idx1).pos;
            const auto e2Pos = data.entry_at(idx2).pos;
            if (reaction->getEducts()[0] == entry1.type) {
                entry1.pos += reaction->getWeight1() * (e2Pos - e1Pos);
            } else {
                entry1.pos += reaction->getWeight1() * (e1Pos - e2Pos);
            }
            fixPos(entry1.pos);
            entry1.type = reaction->getProducts()[0];
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
#endif //READDY_MAIN_SCPUREACTIONUTILS_H
