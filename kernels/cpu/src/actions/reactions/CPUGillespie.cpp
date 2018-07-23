/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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