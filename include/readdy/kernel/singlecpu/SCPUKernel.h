/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
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
