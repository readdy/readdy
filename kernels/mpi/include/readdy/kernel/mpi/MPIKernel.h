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
 * << detailed description >>
 *
 * @file MPIKernel.h
 * @brief Header file for readdy kernel that parallelizes using the message passing interface (MPI)
 * @author chrisfroe
 * @date 28.05.2019
 */

#pragma once


#include <readdy/model/Kernel.h>
#include <readdy/kernel/mpi/MPIStateModel.h>
#include <readdy/kernel/mpi/actions/MPIActionFactory.h>
#include <readdy/kernel/mpi/observables/MPIObservableFactory.h>
#include <readdy/kernel/mpi/model/MPIDomain.h>
#include <readdy/kernel/mpi/Timer.h>

namespace readdy::kernel::mpi {

class MPIKernel : public readdy::model::Kernel {
public:
    static const std::string name;

    MPIKernel();

    ~MPIKernel() override = default;

    explicit MPIKernel(readdy::model::Context ctx);

    MPIKernel(const MPIKernel &) = delete;

    MPIKernel &operator=(const MPIKernel &) = delete;

    MPIKernel(MPIKernel &&) = delete;

    MPIKernel &operator=(MPIKernel &&) = delete;

    // factory method
    static readdy::model::Kernel *create();

    const MPIStateModel &getMPIKernelStateModel() const {
        return _stateModel;
    };

    MPIStateModel &getMPIKernelStateModel() {
        return _stateModel;
    };

    const readdy::model::StateModel &stateModel() const override {
        return _stateModel;
    };

    readdy::model::StateModel &stateModel() override {
        return _stateModel;
    };

    const readdy::model::actions::ActionFactory &actions() const override {
        return _actions;
    };

    readdy::model::actions::ActionFactory &actions() override {
        return _actions;
    };

    const readdy::model::observables::ObservableFactory &observe() const override {
        return _observables;
    };

    readdy::model::observables::ObservableFactory &observe() override {
        return _observables;
    };

    /**
     * Set up domain decomposition.
     * All context and actions must be configured such that cutoffs are known.
     * Consider using SimulationLoop, which always calls initialize(), for simulations.
     */
    void initialize() override;

    const readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() const override {
        return nullptr;
    };

    readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() override {
        return nullptr;
    };

    bool supportsGillespie() const override {
        return false;
    }

    std::shared_ptr<const model::MPIDomain> domain() const {
        return _domain;
    }

protected:
    int rank;
    int worldSize;
    std::string processorName;

    MPIStateModel::Data _data;
    actions::MPIActionFactory _actions;
    observables::MPIObservableFactory _observables;
    MPIStateModel _stateModel;
    std::shared_ptr<model::MPIDomain> _domain{nullptr}; // construction is delayed until initialize()

    // The communicator for the subgroup of actually used workers, can point to _commIfNotWorld or MPI_COMM_WORLD
    MPI_Comm commUsedRanks = MPI_COMM_WORLD;
};

}

extern "C" const char *name();

extern "C" readdy::model::Kernel *createKernel();

