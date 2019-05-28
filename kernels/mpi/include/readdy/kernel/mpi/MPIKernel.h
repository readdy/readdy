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

#include <readdy/model/Kernel.h>

namespace readdy::kernel::mpi {


class MPIKernel : public readdy::model::Kernel {
public:
    static const std::string name;

    MPIKernel();

    ~MPIKernel() override = default;

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

    const model::StateModel &stateModel() const override {
        return _stateModel;
    };

    model::StateModel &stateModel() override {
        return _stateModel;
    };

    const model::actions::ActionFactory &actions() const override {
        return _actions;
    };

    model::actions::ActionFactory &actions() override {
        return _actions;
    };

    const model::observables::ObservableFactory &observe() const override {
        return _observables;
    };

    model::observables::ObservableFactory &observe() override {
        return _observables;
    };

    void initialize() override;

protected:
    int rank;
    std::string processorName;


    MPIStateModel::data_type _data;
    actions::MPIActionFactory _actions;
    observables::MPIObservableFactory _observables;
    MPIStateModel _stateModel;
};

}

extern "C" const char *name();

extern "C" readdy::model::Kernel *createKernel();

