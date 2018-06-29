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
 * @date 20.10.16
 */
#include <readdy/kernel/cpu/actions/reactions/CPUGillespie.h>


namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

CPUGillespie::CPUGillespie(CPUKernel *const kernel, scalar timeStep) : super(timeStep), kernel(kernel) {}

void CPUGillespie::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    const auto &ctx = kernel->context();
    if(ctx.reactions().nOrder1() == 0 && ctx.reactions().nOrder2() == 0) {
        return;
    }
    auto &stateModel = kernel->getCPUKernelStateModel();
    auto data = stateModel.getParticleData();
    const auto nl = stateModel.getNeighborList();

    if(ctx.recordReactionCounts()) {
        stateModel.resetReactionCounts();
    }

    scalar alpha = 0.0;
    std::vector<event_t> events;
    gatherEvents(kernel, readdy::util::range<event_t::index_type>(0, data->size()), nl, data, alpha, events);
    if(ctx.recordReactionsWithPositions()) {
        stateModel.reactionRecords().clear();
        if(ctx.recordReactionCounts()) {
            auto &counts = stateModel.reactionCounts();
            auto particlesUpdate = handleEventsGillespie(kernel, timeStep(), false, false, std::move(events),
                                                         &stateModel.reactionRecords(), &counts);
            data->update(std::move(particlesUpdate));
        } else {
            auto particlesUpdate = handleEventsGillespie(kernel, timeStep(), false, false, std::move(events),
                                                         &stateModel.reactionRecords(), nullptr);
            data->update(std::move(particlesUpdate));
        }
    } else {
        if(ctx.recordReactionCounts()) {
            auto &counts = stateModel.reactionCounts();
            auto particlesUpdate = handleEventsGillespie(kernel, timeStep(), false, false, std::move(events),
                                                         nullptr, &counts);
            data->update(std::move(particlesUpdate));
        } else {
            auto particlesUpdate = handleEventsGillespie(kernel, timeStep(), false, false, std::move(events),
                                                         nullptr, nullptr);
            data->update(std::move(particlesUpdate));
        }
    }
}

}
}
}
}
}