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
 * @brief Implementation of CPU program Compartments
 * @author chrisfroe
 * @date 18.10.16
 */

#include <readdy/kernel/cpu/actions/CPUEvaluateCompartments.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {

CPUEvaluateCompartments::CPUEvaluateCompartments(CPUKernel *const kernel) : kernel(kernel) {}

void CPUEvaluateCompartments::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    const auto &ctx = kernel->context();
    const auto &compartments = ctx.compartments().get();
    for(auto& e : *kernel->getCPUKernelStateModel().getParticleData()) {
        if(!e.deactivated) {
            for (const auto &compartment : compartments) {
                if (compartment->isContained(e.pos)) {
                    const auto &conversions = compartment->getConversions();
                    const auto convIt = conversions.find(e.type);
                    if (convIt != conversions.end()) {
                        e.type = (*convIt).second;
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
