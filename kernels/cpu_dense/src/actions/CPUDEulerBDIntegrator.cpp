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
 * @file EulerBDIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/common/thread/scoped_thread.h>
#include <readdy/kernel/cpu_dense/actions/CPUDEulerBDIntegrator.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace actions {

namespace rnd = readdy::model::rnd;
namespace thd = readdy::util::thread;

void CPUDEulerBDIntegrator::perform() {
    auto& pd = *kernel->getCPUDKernelStateModel().getParticleData();
    const auto size = pd.size();
    std::vector<thd::scoped_thread> threads;
    threads.reserve(kernel->getNThreads());
    const std::size_t grainSize = size / kernel->getNThreads();

    const auto &context = kernel->getKernelContext();
    using iter_t = decltype(pd.begin());

    const auto dt = timeStep;

    auto worker = [&context, &pd, dt](iter_t entry_begin, iter_t entry_end)  {
        const auto &fixPos = context.getFixPositionFun();
        const auto kbt = context.getKBT();
        for (iter_t it = entry_begin; it != entry_end; ++it) {
            const double D = context.particleTypeRegistry().getDiffusionConstant(it->type);
            const auto randomDisplacement = std::sqrt(2. * D * dt) * rnd::normal3(0, 1);
            const auto deterministicDisplacement = it->force * dt * D / kbt;
            pd.displace(*it, randomDisplacement + deterministicDisplacement);
        }

    };
    auto work_iter = pd.begin();
    {
        for (unsigned int i = 0; i < kernel->getNThreads() - 1; ++i) {
            threads.push_back(thd::scoped_thread(std::thread(worker, work_iter, work_iter + grainSize)));
            work_iter += grainSize;
        }
        threads.push_back(thd::scoped_thread(std::thread(worker, work_iter, pd.end())));
    }

}

CPUDEulerBDIntegrator::CPUDEulerBDIntegrator(CPUDKernel *const kernel, double timeStep)
        : super(timeStep), kernel(kernel) {}

}
}
}
}