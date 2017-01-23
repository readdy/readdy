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
 * @file Gillespie.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */

#include <readdy/kernel/cpu_dense/actions/reactions/CPUDGillespie.h>


namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace actions {
namespace reactions {

CPUDGillespie::CPUDGillespie(const CPUDKernel *const kernel, double timeStep) : super(timeStep), kernel(kernel) {}

void CPUDGillespie::perform() {
    const auto &ctx = kernel->getKernelContext();
    auto data = kernel->getKernelStateModel().getParticleData();
    const auto &dist = ctx.getDistSquaredFun();
    const auto &fixPos = ctx.getFixPositionFun();
    const auto nl = kernel->getKernelStateModel().getNeighborList();

    double alpha = 0.0;
    std::vector<event_t> events;
    gatherEvents(kernel, readdy::util::range<event_t::index_type>(0, data->size()),
                 nl, *data, alpha, events, dist);
    auto particlesUpdate = handleEventsGillespie(kernel, timeStep, false, true, std::move(events));

    // update data structure
    data->deactivateMarked();
    data->update(std::move(particlesUpdate));
}

}
}
}
}
}