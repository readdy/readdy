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
 * The SingleCPUCompartments program defines compartments via characteristic functions that map from Vec3 to bool. For every compartment
 * one can then define conversions that should take place as soon as a particle enters the compartment. Note that the user is responsible for
 * keeping the compartments disjoint.
 *
 * @file SingleCPUCompartments.cpp
 * @brief SingleCPUCompartments defines compartments via characteristic functions and performs conversions of particle species on these compartments.
 * @author chrisfroe
 * @date 13.10.16
 */
#pragma once
#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {

class SCPUEvaluateCompartments : public readdy::model::actions::EvaluateCompartments {
public:
    explicit SCPUEvaluateCompartments(SCPUKernel* kernel);

    void perform(const util::PerformanceNode &node) override;

protected:
    SCPUKernel *const kernel;
};

}
}
}
}
