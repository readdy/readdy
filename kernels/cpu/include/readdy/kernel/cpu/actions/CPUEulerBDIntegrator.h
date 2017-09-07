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
 * @file CPUEulerBDIntegrator.h
 * @brief << brief description >>
 * @author clonker
 * @date 07.07.16
 */

#pragma once
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/model/actions/Actions.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
class CPUEulerBDIntegrator : public readdy::model::actions::EulerBDIntegrator {

public:
    CPUEulerBDIntegrator(CPUKernel *kernel, readdy::scalar timeStep);

    void perform(const util::PerformanceNode &node) override;

private:
    CPUKernel *kernel;
};
}
}
}
}
