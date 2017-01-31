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
 * @file Compartments.cpp
 * @brief Implementation of SingleCPU program Compartments
 * @author chrisfroe
 * @date 13.10.16
 */

#include <readdy/kernel/singlecpu/actions/SCPUEvaluateCompartments.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {

SCPUEvaluateCompartments::SCPUEvaluateCompartments(SCPUKernel *const kernel) : kernel(kernel) {}

void SCPUEvaluateCompartments::perform() {
    const auto &ctx = kernel->getKernelContext();
    const auto & compartments = ctx.getCompartments();
    auto data = kernel->getSCPUKernelStateModel().getParticleData();
    const auto &compartments = ctx.getCompartments();
    for(auto& entry : *data) {
        if(!entry.is_deactivated()) {
            for (auto i=0; i<compartments.size(); ++i) {
                if (compartments[i]->isContained(entry.position())) {
                    const auto &conversions = compartments[i]->getConversions();
                    const auto convIt = conversions.find(entry.type);
                    if (convIt != conversions.end()) {
                        entry.type = convIt->second;
                    }
                }
            }
        }
    }
}

}
}
}
}

