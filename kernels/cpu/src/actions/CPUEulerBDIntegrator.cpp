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
 * @file CPUEulerBDIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.07.16
 */

#include <readdy/kernel/cpu/actions/CPUEulerBDIntegrator.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {

namespace rnd = readdy::model::rnd;

void CPUEulerBDIntegrator::perform(const readdy::util::PerformanceNode &node) {
    auto t = node.timeit();
    auto data = kernel->getCPUKernelStateModel().getParticleData();
    const auto size = data->size();

    const auto &context = kernel->context();
    using iter_t = data::EntryDataContainer::iterator;

    const auto dt = timeStep();

    auto worker = [&context, data, dt](std::size_t, std::size_t beginIdx, iter_t entry_begin, iter_t entry_end)  {
        const auto kbt = context.kBT();
        std::size_t idx = beginIdx;
        const auto &box = context.boxSize().data();
        const auto &pbc = context.periodicBoundaryConditions().data();
        for (auto it = entry_begin; it != entry_end; ++it, ++idx) {
            if(!it->deactivated) {
                const scalar D = context.particleTypes().diffusionConstantOf(it->type);
                const auto randomDisplacement = std::sqrt(2. * D * dt) * rnd::normal3<readdy::scalar>(0, 1);
                const auto deterministicDisplacement = it->force * dt * D / kbt;
                it->pos += randomDisplacement + deterministicDisplacement;
                bcs::fixPosition(it->pos, box, pbc);
            }
        }
    };

    std::vector<util::thread::joining_future<void>> waitingFutures;
    waitingFutures.reserve(kernel->getNThreads());
    auto &pool  = kernel->pool();
    {
        auto it = data->begin();
        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(kernel->getNThreads());

        auto granularity = kernel->getNThreads();
        const std::size_t grainSize = size / granularity;

        std::size_t idx = 0;
        for (auto i = 0_z; i < granularity-1; ++i) {
            auto itNext = it + grainSize;
            if(it != itNext) {
                waitingFutures.emplace_back(pool.push(worker, idx, it, itNext));
            }
            it = itNext;
            idx += grainSize;
        }
        if(it != data->end()) {
            waitingFutures.emplace_back(pool.push(worker, idx, it, data->end()));
        }
    }

}

CPUEulerBDIntegrator::CPUEulerBDIntegrator(CPUKernel *kernel, scalar timeStep)
        : readdy::model::actions::EulerBDIntegrator(timeStep), kernel(kernel) {}

}
}
}
}
