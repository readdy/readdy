/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
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
 * « detailed description »
 *
 * @file MPIActions.h
 * @brief « brief description »
 * @author chrisfroe
 * @date 03.06.19
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/mpi/MPIKernel.h>

namespace readdy::kernel::mpi::actions {

class MPIEulerBDIntegrator : public readdy::model::actions::EulerBDIntegrator {
public:
    MPIEulerBDIntegrator(MPIKernel *kernel, readdy::scalar timeStep) : kernel(kernel), EulerBDIntegrator(timeStep) {}

    void perform() override;

private:
    MPIKernel *kernel;
};

class MPICalculateForces : public readdy::model::actions::CalculateForces {
public:
    explicit MPICalculateForces(MPIKernel *kernel) : CalculateForces(), kernel(kernel) {}

    void perform() override {
        const auto &context = kernel->context();
        if (context.recordVirial()) {
            performImpl<true>();
        } else {
            performImpl<false>();
        }
    }

private:
    MPIKernel *kernel;

    template<bool COMPUTE_VIRIAL>
    void performImpl();
};

class MPICreateNeighborList : public readdy::model::actions::CreateNeighborList {
public:
    MPICreateNeighborList(MPIKernel *kernel, scalar cutoffDistance) : CreateNeighborList(cutoffDistance),
                                                                      kernel(kernel) {}

    void perform() override {
        kernel->getMPIKernelStateModel().initializeNeighborList(cutoffDistance());
    }

private:
    MPIKernel *kernel;

};

class MPIUpdateNeighborList : public readdy::model::actions::UpdateNeighborList {
public:
    explicit MPIUpdateNeighborList(MPIKernel *kernel) : UpdateNeighborList(), kernel(kernel) {}

    void perform() override {
        kernel->getMPIKernelStateModel().updateNeighborList();
    }

private:
    MPIKernel *kernel;
};

class MPIClearNeighborList : public readdy::model::actions::ClearNeighborList {
public:
    explicit MPIClearNeighborList(MPIKernel *kernel) : ClearNeighborList(), kernel(kernel) {}

    void perform() override {
        kernel->getMPIKernelStateModel().clearNeighborList();
    }

private:
    MPIKernel *kernel;
};

class MPIEvaluateCompartments : public readdy::model::actions::EvaluateCompartments {
public:
    explicit MPIEvaluateCompartments(MPIKernel *kernel) : EvaluateCompartments(), kernel(kernel) {}

    void perform() override {
        const auto &ctx = kernel->context();
        const auto &compartments = ctx.compartments().get();
        for (auto &e : *kernel->getMPIKernelStateModel().getParticleData()) {
            if (!e.deactivated) {
                for (const auto &compartment : compartments) {
                    if (compartment->isContained(e.pos)) {
                        const auto &conversions = compartment->getConversions();
                        const auto convIt = conversions.find(e.type);
                        if (convIt != conversions.end()) {
                            e.type = (*convIt).second;
                        }
                    }
                }
            }
        }
    }

protected:
    MPIKernel *kernel;
};

namespace reactions {

class MPIUncontrolledApproximation : public readdy::model::actions::reactions::UncontrolledApproximation {
public:
    MPIUncontrolledApproximation(MPIKernel *kernel, readdy::scalar timeStep) : UncontrolledApproximation(timeStep),
                                                                               kernel(kernel) {}

    void perform() override;

protected:
    MPIKernel *const kernel;
};
}

}