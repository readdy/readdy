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
 * @file IOUtils.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/IOUtils.h>

namespace readdy {
namespace model {
void writeReactionInformation(io::Group &group, const KernelContext &context) {
    auto subgroup = group.createGroup("./registered_reactions");
    auto groupO1 = subgroup.createGroup("./order1");
    auto groupO2 = subgroup.createGroup("./order2");
    for (const auto &t1 : context.getAllRegisteredParticleTypes()) {
        const auto &rO1 = context.getOrder1Reactions(t1);
        std::vector<std::string> lO1;
        lO1.reserve(rO1.size());
        std::for_each(rO1.begin(), rO1.end(), [&lO1](reactions::Reaction<1> *r) { lO1.push_back(r->getName()); });
        if (!lO1.empty()) {
            groupO1.write(context.getParticleName(t1) + "[id=" + std::to_string(t1) + "]", lO1);
        }
        for (const auto &t2 : context.getAllRegisteredParticleTypes()) {
            if (t2 < t1) continue;
            const auto &rO2 = context.getOrder2Reactions(t1, t2);
            std::vector<std::string> lO2;
            lO2.reserve(rO2.size());
            std::for_each(rO2.begin(), rO2.end(), [&lO2](reactions::Reaction<2> *r) { lO2.push_back(r->getName()); });
            if (!lO2.empty()) {
                groupO2.write(
                        context.getParticleName(t1) + "[id=" + std::to_string(t1) + "] + " +
                        context.getParticleName(t2) + "[id=" + std::to_string(t2) + "]", lO2);
            }
        }
    }
}

}
}