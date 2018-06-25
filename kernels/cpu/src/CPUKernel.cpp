/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * @file CPUKernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12/11/17
 */

#include <readdy/kernel/cpu/CPUKernel.h>


namespace readdy {
namespace kernel {
namespace cpu {

const std::string CPUKernel::name = "CPU";

readdy::model::Kernel *CPUKernel::create() {
    return new CPUKernel();
}

CPUKernel::CPUKernel() : readdy::model::Kernel(name), _pool(readdy_default_n_threads()),
                         _data(_context, _pool), _actions(this),
                         _observables(this), _topologyActionFactory(_context, _data),
                         _stateModel(_data, _context, _pool, &_topologyActionFactory) {}

void CPUKernel::initialize() {
    readdy::model::Kernel::initialize();

    const auto &fullConfiguration = context().kernelConfiguration();

    const auto &configuration = fullConfiguration.cpu;
    // thread config
    setNThreads(static_cast<std::uint32_t>(configuration.threadConfig.getNThreads()));
    {
        // state model config
        _stateModel.configure(configuration);
    }
    for (auto &top : _stateModel.topologies()) {
        top->configure();
        top->updateReactionRates(context().topologyRegistry().structuralReactionsOf(top->type()));
    }
    _stateModel.reactionRecords().clear();
    _stateModel.resetReactionCounts();
    _stateModel.virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};
}

}
}
}

const char *name() {
    return readdy::kernel::cpu::CPUKernel::name.c_str();
}

readdy::model::Kernel *createKernel() {
    return readdy::kernel::cpu::CPUKernel::create();
}
