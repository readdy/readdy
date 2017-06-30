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

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {

namespace rnd = readdy::model::rnd;
namespace thd = readdy::util::thread;

void CPUEulerBDIntegrator::perform() {
    auto& pd = *kernel->getCPUKernelStateModel().getParticleData();
    const auto size = pd.size();

    const auto &context = kernel->getKernelContext();
    using iter_t = decltype(pd.begin());

    const auto dt = timeStep;

    auto worker = [&context, &pd, dt](std::size_t id, iter_t entry_begin, iter_t entry_end)  {
        const auto &fixPos = context.getFixPositionFun();
        const auto kbt = context.getKBT();
        for (iter_t it = entry_begin; it != entry_end; ++it) {
            if(!it->is_deactivated()) {
                const double D = context.particle_types().diffusion_constant_of(it->type);
                const auto randomDisplacement = std::sqrt(2. * D * dt) * rnd::normal3(0, 1);
                const auto deterministicDisplacement = it->force * dt * D / kbt;
                pd.displace(*it, randomDisplacement + deterministicDisplacement);
            }
        }
    };

    auto work_iter = pd.begin();
    {
        auto& executor = kernel->executor();
        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(kernel->getNThreads());

        auto granularity = kernel->getNThreads();
        const std::size_t grainSize = size / granularity;

        for (unsigned int i = 0; i < granularity - 1; ++i) {
            executables.push_back(executor.pack(worker, work_iter, work_iter+grainSize));
            work_iter += grainSize;
        }
        executables.push_back(executor.pack(worker, work_iter, pd.end()));
        executor.execute_and_wait(std::move(executables));
    }

}

CPUEulerBDIntegrator::CPUEulerBDIntegrator(CPUKernel *kernel, double timeStep)
        : readdy::model::actions::EulerBDIntegrator(timeStep), kernel(kernel) {}

}
}
}
}