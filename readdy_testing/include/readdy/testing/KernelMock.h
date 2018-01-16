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


/**
 * << detailed description >>
 *
 * @file KernelMock.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_MAIN_KERNELMOCK_H
#define READDY_MAIN_KERNELMOCK_H

#include <readdy/model/Kernel.h>
#include <gmock/gmock.h>

namespace readdy {
namespace testing {

class FakeActionFactory : public readdy::model::actions::ActionFactory {
public:
    std::unique_ptr<model::actions::AddParticles>
    addParticles(const std::vector<model::Particle> &particles) const override {
        return nullptr;
    }

    std::unique_ptr<model::actions::EulerBDIntegrator> eulerBDIntegrator(scalar timeStep) const override {
        return nullptr;
    }

    std::unique_ptr<model::actions::CalculateForces> calculateForces() const override {
        return nullptr;
    }

    std::unique_ptr<model::actions::UpdateNeighborList>
    updateNeighborList(model::actions::UpdateNeighborList::Operation operation, scalar skinSize) const override {
        return nullptr;
    }

    std::unique_ptr<model::actions::EvaluateCompartments> evaluateCompartments() const override {
        return nullptr;
    }

    std::unique_ptr<model::actions::reactions::UncontrolledApproximation>
    uncontrolledApproximation(scalar timeStep) const override {
        return nullptr;
    }

    std::unique_ptr<model::actions::reactions::Gillespie> gillespie(scalar timeStep) const override {
        return nullptr;
    }

    std::unique_ptr<model::actions::top::EvaluateTopologyReactions>
    evaluateTopologyReactions(scalar timeStep) const override {
        return nullptr;
    }
};

class KernelMock : public readdy::model::Kernel {

public:
    explicit KernelMock(const std::string &name) : Kernel(name) {}

    MOCK_METHOD0(actions, readdy::model::actions::ActionFactory & (void));

    MOCK_CONST_METHOD0(actions, const readdy::model::actions::ActionFactory & (
            void));

    MOCK_METHOD0(stateModel, readdy::model::StateModel &(void));
    MOCK_CONST_METHOD0(stateModel, const readdy::model::StateModel & (void));

    MOCK_METHOD0(observe, readdy::model::observables::ObservableFactory &(void));
    MOCK_CONST_METHOD0(observe, const readdy::model::observables::ObservableFactory & (void));

    MOCK_METHOD0(getTopologyActionFactory, readdy::model::top::TopologyActionFactory*const (void));
    MOCK_CONST_METHOD0(getTopologyActionFactory, const readdy::model::top::TopologyActionFactory*const (void));
};
}
}
#endif //READDY_MAIN_KERNELMOCK_H
