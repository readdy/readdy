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
    for (const auto &t1 : context.particle_types().types_flat()) {
        const auto &reactionsOrder1 = context.reactions().order1_by_type(t1);
        std::vector<std::string> labelsOrder1;
        labelsOrder1.reserve(reactionsOrder1.size());
        std::for_each(reactionsOrder1.begin(), reactionsOrder1.end(),
                      [&labelsOrder1](reactions::Reaction<1> *r) { labelsOrder1.push_back(r->getName()); });
        if (!labelsOrder1.empty()) {
            groupO1.write(context.particle_types().name_of(t1) + "[id=" + std::to_string(t1) + "]", labelsOrder1);
        }
        for (const auto &t2 : context.particle_types().types_flat()) {
            if (t2 < t1) continue;
            const auto &reactionsOrder2 = context.reactions().order2_by_type(t1, t2);
            std::vector<std::string> labelsOrder2;
            labelsOrder2.reserve(reactionsOrder2.size());
            std::for_each(reactionsOrder2.begin(), reactionsOrder2.end(),
                          [&labelsOrder2](reactions::Reaction<2> *r) { labelsOrder2.push_back(r->getName()); });
            if (!labelsOrder2.empty()) {
                groupO2.write(
                        context.particle_types().name_of(t1) + "[id=" + std::to_string(t1) + "] + " +
                        context.particle_types().name_of(t2) + "[id=" + std::to_string(t2) + "]", labelsOrder2);
            }
        }
    }
}

void writeParticleTypeInformation(io::Group &group, const KernelContext &context) {
    auto subgroup = group.createGroup("./particle_types");
    const auto& types = context.particle_types().type_mapping();
    for(auto it = types.begin(); it != types.end(); ++it) {
        std::vector<std::string> strvec(1);
        strvec.at(0) = it->first;
        subgroup.write(std::to_string(it->second), it->first);
    }
}

}
}