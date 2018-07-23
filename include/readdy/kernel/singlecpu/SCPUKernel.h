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

#pragma once

#include <readdy/model/RandomProvider.h>
#include <readdy/model/Kernel.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservableFactory.h>
#include <readdy/kernel/singlecpu/model/topologies/SCPUTopologyActionFactory.h>
#include <readdy/kernel/singlecpu/actions/SCPUActionFactory.h>

namespace readdy {
namespace kernel {
namespace scpu {

class SCPUKernel : public readdy::model::Kernel {
public:

    static const std::string name;

    SCPUKernel();

    ~SCPUKernel() override;

    // factory method
    static std::unique_ptr<SCPUKernel> create();

    const SCPUStateModel &getSCPUKernelStateModel() const {
        return _model;
    };

    SCPUStateModel &getSCPUKernelStateModel() {
        return _model;
    };

    void initialize() override;

    const readdy::model::observables::ObservableFactory &observe() const override {
        return _observables;
    };

    readdy::model::observables::ObservableFactory &observe() override {
        return _observables;
    };

    const readdy::model::actions::ActionFactory &actions() const override {
        return _actionFactory;
    };

    readdy::model::actions::ActionFactory &actions() override {
        return _actionFactory;
    };

    const readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() const override {
        return &_topologyActionFactory;
    };

    readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() override {
        return &_topologyActionFactory;
    };

    const readdy::model::StateModel &stateModel() const override {
        return _model;
    };

    readdy::model::StateModel &stateModel() override {
        return _model;
    };

private:
    model::top::SCPUTopologyActionFactory _topologyActionFactory;
    SCPUStateModel _model;
    actions::SCPUActionFactory _actionFactory;
    observables::SCPUObservableFactory _observables;
};

}
}
}
