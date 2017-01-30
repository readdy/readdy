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

CPUEvaluateCompartments::CPUEvaluateCompartments(const CPUKernel *const kernel) : kernel(kernel) {}

void CPUEvaluateCompartments::perform() {
    const auto &ctx = kernel->getKernelContext();
    const auto &compartments = ctx.getCompartments();
    for(auto& e : *kernel->getKernelStateModel().getParticleData()) {
        if(!e.is_deactivated()) {
            for (auto i = 0; i < compartments.size(); ++i) {
                if (compartments[i]->isContained(e.position())) {
                    const auto conversions = compartments[i]->getConversions();
                    if (conversions.find(e.type) != conversions.end()) {
                        e.type = conversions.at(e.type);
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
