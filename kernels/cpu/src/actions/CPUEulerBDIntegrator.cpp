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
 * @file CPUEulerBDIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.07.16
 */

#include <readdy/kernel/cpu/actions/CPUEulerBDIntegrator.h>
#include <readdy/kernel/cpu/data/NLDataContainer.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {

namespace rnd = readdy::model::rnd;

void CPUEulerBDIntegrator::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    auto data = kernel->getCPUKernelStateModel().getParticleData();
    const auto size = data->size();

    const auto &context = kernel->context();
    using iter_t = data::EntryDataContainer::iterator;

    const auto dt = timeStep;

    auto worker = [&context, data, dt](std::size_t id, std::size_t beginIdx,  iter_t entry_begin, iter_t entry_end)  {
        const auto &fixPos = context.fixPositionFun();
        const auto kbt = context.kBT();
        std::size_t idx = beginIdx;
        for (auto it = entry_begin; it != entry_end; ++it, ++idx) {
            if(!it->deactivated) {
                const scalar D = context.particle_types().diffusionConstantOf(it->type);
                const auto randomDisplacement = std::sqrt(2. * D * dt) * rnd::normal3<readdy::scalar>(0, 1);
                const auto deterministicDisplacement = it->force * dt * D / kbt;
                data->displace(idx, randomDisplacement + deterministicDisplacement);
            }
        }
    };

    auto work_iter = data->begin();
    {
        auto& executor = kernel->executor();
        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(kernel->getNThreads());

        auto granularity = kernel->getNThreads();
        const std::size_t grainSize = size / granularity;

        std::size_t idx = 0;
        for (unsigned int i = 0; i < granularity - 1; ++i) {
            executables.push_back(executor.pack(worker, idx, work_iter, work_iter+grainSize));
            work_iter += grainSize;
            idx += grainSize;
        }
        executables.push_back(executor.pack(worker, idx, work_iter, data->end()));
        executor.execute_and_wait(std::move(executables));
    }

}

CPUEulerBDIntegrator::CPUEulerBDIntegrator(CPUKernel *kernel, scalar timeStep)
        : readdy::model::actions::EulerBDIntegrator(timeStep), kernel(kernel) {}

}
}
}
}