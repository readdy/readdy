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
 * @file SCPUCreateNeighborList.cpp
 * @brief Implementation of CreateNeighborList action for single CPU kernel
 * @author clonker
 * @author chrisfroe
 * @date 11.07.16
 */

#include <readdy/kernel/singlecpu/actions/SCPUCreateNeighborList.h>

namespace readdy::kernel::scpu::actions {

void SCPUCreateNeighborList::perform() {
    kernel->getSCPUKernelStateModel().initializeNeighborList(_cutoffDistance);
}

SCPUCreateNeighborList::SCPUCreateNeighborList(SCPUKernel *kernel, scalar cutoffDistance)
        : CreateNeighborList(cutoffDistance), kernel(kernel) {
}

SCPUUpdateNeighborList::SCPUUpdateNeighborList(SCPUKernel *kernel) : kernel(kernel) {}

void SCPUUpdateNeighborList::perform() {
    kernel->getSCPUKernelStateModel().updateNeighborList();
}

SCPUClearNeighborList::SCPUClearNeighborList(SCPUKernel *kernel) : kernel(kernel) {}

void SCPUClearNeighborList::perform() {
    kernel->getSCPUKernelStateModel().clearNeighborList();
}

}
