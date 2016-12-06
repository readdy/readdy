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
 * @file SingleCPUReactionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/kernel/singlecpu/reactions/SCPUReactionFactory.h>
#include <readdy/kernel/singlecpu/reactions/SCPUReactions.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace reactions {
SCPUReactionFactory::SCPUReactionFactory(SCPUKernel const *const kernel) : kernel(
        kernel) {

}

readdy::model::reactions::Conversion *
SCPUReactionFactory::createConversion(const std::string &name, unsigned int from, unsigned int to,
                                           const double rate) const {
    return new SCPUConversion(name, from, to, rate);
}

readdy::model::reactions::Enzymatic *
SCPUReactionFactory::createEnzymatic(const std::string &name, unsigned int catalyst,
                                          unsigned int from, unsigned int to, const double rate,
                                          const double eductDistance) const {
    return new SCPUEnzymatic(name, catalyst, from, to, rate, eductDistance);
}

readdy::model::reactions::Fission *
SCPUReactionFactory::createFission(const std::string &name, unsigned int from, unsigned int to1,
                                        unsigned int to2,
                                        const double rate, const double productDistance,
                                        const double weight1, const double weight2) const {
    return new SCPUFission(name, from, to1, to2, rate, productDistance, weight1, weight2);
}

readdy::model::reactions::Fusion *
SCPUReactionFactory::createFusion(const std::string &name, unsigned int from1, unsigned int from2,
                                       unsigned int to, const double rate,
                                       const double eductDistance, const double weight1,
                                       const double weight2) const {
    return new SCPUFusion(name, from1, from2, to, rate, eductDistance, weight1, weight2);
}

}
}
}
}