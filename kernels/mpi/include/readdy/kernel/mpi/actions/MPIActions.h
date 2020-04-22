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

#include <utility>
#include <readdy/api/Saver.h>

namespace readdy::kernel::mpi::actions {

class MPIAddParticles : public readdy::model::actions::AddParticles {
public:
    MPIAddParticles(MPIKernel *kernel, const std::vector<readdy::model::Particle> &particles) : AddParticles(kernel, particles), kernel(kernel) {}

    void perform() override {
        if (kernel) {
            kernel->getMPIKernelStateModel().distributeParticles(particles);
        } else {
            throw std::runtime_error("invalid kernel");
        }
    }
private:
    MPIKernel *const kernel;
};

class MPIEulerBDIntegrator : public readdy::model::actions::EulerBDIntegrator {
public:
    MPIEulerBDIntegrator(MPIKernel *kernel, readdy::scalar timeStep) : kernel(kernel), EulerBDIntegrator(timeStep) {}

    void perform() override;

private:
    MPIKernel *const kernel;
};

class MPICalculateForces : public readdy::model::actions::CalculateForces {
public:
    explicit MPICalculateForces(MPIKernel *kernel) : CalculateForces(), kernel(kernel) {}

    void perform() override {
        if (kernel->domain().isWorkerRank()) {
            const auto &context = kernel->context();
            if (context.recordVirial()) {
                performImpl<true>();
            } else {
                performImpl<false>();
            }
        } else {
            readdy::log::trace("MPICalculateForces::perform is noop for non workers");
        }
    }

private:
    MPIKernel *const kernel;

    template<bool COMPUTE_VIRIAL>
    void performImpl();
};

/** no-op but needed to run the default simulation loop */
class MPICreateNeighborList : public readdy::model::actions::CreateNeighborList {
public:
    MPICreateNeighborList(MPIKernel *kernel) : CreateNeighborList(kernel->context().calculateMaxCutoff()) {}
    void perform() override {/* no-op, this neighborlist is initialized by construction */}
};

class MPIUpdateNeighborList : public readdy::model::actions::UpdateNeighborList {
public:
    explicit MPIUpdateNeighborList(MPIKernel *kernel) : UpdateNeighborList(), kernel(kernel) {}

    void perform() override {
        if (kernel->domain().isWorkerRank()) {
            kernel->getMPIKernelStateModel().updateNeighborList();
        } else {
            readdy::log::trace("MPIUpdateNeighborList::perform is noop for non workers");
        }
    }

private:
    MPIKernel *const kernel;
};

/** no-op but needed to run the default simulation loop */
class MPIClearNeighborList : public readdy::model::actions::ClearNeighborList {
public:
    explicit MPIClearNeighborList() : ClearNeighborList() {}
    void perform() override {/* no-op */}
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
    MPIKernel *const kernel;
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

class MPIEvaluateObservables : public readdy::model::actions::EvaluateObservables {
public:
    explicit MPIEvaluateObservables(MPIKernel *kernel) : kernel(kernel) {}

    void perform(TimeStep t) override {
        kernel->evaluateObservables(t);
    }

private:
    MPIKernel *kernel;
};

class MPIMakeCheckpoint : public readdy::model::actions::MakeCheckpoint {
public:
    MPIMakeCheckpoint(MPIKernel *kernel, const std::string& base, std::size_t maxNSaves,
                      const std::string &checkpointFormat)
            : kernel(kernel), saver(base, maxNSaves, checkpointFormat) {}

    void perform(TimeStep t) override {
        // todo sync (MPIGather) the state to master's stateModel, then makeCheckpoint as usual and clear stateModel
        if (kernel->domain().isMasterRank()) {
            // todo gather
            // make checkpoint on master
            //saver.makeCheckpoint(kernel, t);
            // clear state model
            //kernel->getMPIKernelStateModel().clear();
        } else if (kernel->domain().isWorkerRank()) {
            // todo gather, send responsible particles
        } else {
            // no op for idlers
        }
    }

    std::string describe() const override {
        return saver.describe();
    }
private:
    MPIKernel *kernel;
    readdy::api::Saver saver;
};

class MPIInitializeKernel : public readdy::model::actions::InitializeKernel {
public:
    MPIInitializeKernel(MPIKernel *kernel) : kernel(kernel) {}

    void perform() override {
        kernel->initialize();
    }
private:
    MPIKernel *kernel;
};

}
