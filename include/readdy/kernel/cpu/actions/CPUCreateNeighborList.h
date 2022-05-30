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
 * @file CPUNeighborList.h
 * @brief CPU kernel declaration of NeighborList Action
 * @author clonker
 * @date 13.07.16
 */

#pragma once
#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy::kernel::cpu::actions {

class CPUCreateNeighborList : public readdy::model::actions::CreateNeighborList {
public:
    CPUCreateNeighborList(CPUKernel *kernel, scalar cutoffDistance)
            : CreateNeighborList(cutoffDistance), kernel(kernel) {}

    void perform() override {
        kernel->getCPUKernelStateModel().initializeNeighborList(_cutoffDistance);
    }

private:
    CPUKernel *const kernel;
};

class CPUUpdateNeighborList : public readdy::model::actions::UpdateNeighborList {
public:
    explicit CPUUpdateNeighborList(CPUKernel *kernel) : kernel(kernel) {}

    void perform() override {
        kernel->getCPUKernelStateModel().updateNeighborList();
    }

private:
    CPUKernel *const kernel;
};

class CPUClearNeighborList : public readdy::model::actions::ClearNeighborList {
public:
    explicit CPUClearNeighborList(CPUKernel *kernel) : kernel(kernel) {}

    void perform() override {
        kernel->getCPUKernelStateModel().clearNeighborList();
    }

private:
    CPUKernel *const kernel;
};

}
