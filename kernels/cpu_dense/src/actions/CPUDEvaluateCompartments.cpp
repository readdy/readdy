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
 * @file Compartments.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/actions/CPUDEvaluateCompartments.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace actions {

CPUDEvaluateCompartments::CPUDEvaluateCompartments(CPUDKernel *const kernel) : kernel(kernel) {}

void CPUDEvaluateCompartments::perform() {
    const auto &ctx = kernel->getKernelContext();
    const auto &compartments = ctx.getCompartments();
    for(auto& e : *kernel->getKernelStateModel().getParticleData()) {
        for (auto i = 0; i < compartments.size(); ++i) {
            if (compartments[i]->isContained(e.position())) {
                const auto &conversions = compartments[i]->getConversions();
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
