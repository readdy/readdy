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


//
// Created by clonker on 07.03.16.
//

#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
const std::string SCPUKernel::name = "SingleCPU";

SCPUKernel::SCPUKernel() : readdy::model::Kernel(name), _actionFactory(this), _topologyActionFactory(this),
                           _model(_context, &_topologyActionFactory), _observables(this) {}

/**
 * factory method
 */
std::unique_ptr<SCPUKernel> SCPUKernel::create() {
    return std::make_unique<SCPUKernel>();
}

/**
 * Destructor: default
 */
SCPUKernel::~SCPUKernel() = default;

void SCPUKernel::initialize() {
    readdy::model::Kernel::initialize();
    for(auto& top : getSCPUKernelStateModel().topologies()) {
        top->configure();
        top->updateReactionRates(context().topologyRegistry().structuralReactionsOf(top->type()));
    }
    getSCPUKernelStateModel().reactionRecords().clear();
    getSCPUKernelStateModel().resetReactionCounts();
    getSCPUKernelStateModel().virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};
}

}
}
}




