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

void SCPUEvaluateCompartments::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    const auto &ctx = kernel->context();
    const auto & compartments = ctx.compartments().get();
    auto data = kernel->getSCPUKernelStateModel().getParticleData();
    for(auto& entry : *data) {
        if(!entry.is_deactivated()) {
            for (const auto &compartment : compartments) {
                if (compartment->isContained(entry.position())) {
                    const auto& conversions = compartment->getConversions();
                    auto it = conversions.find(entry.type);
                    if (it != conversions.end()) {
                        entry.type = it->second;
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

