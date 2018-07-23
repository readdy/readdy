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
