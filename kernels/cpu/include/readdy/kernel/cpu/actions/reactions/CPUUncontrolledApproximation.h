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
 * @file UncontrolledApproximation.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#pragma once
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

class CPUUncontrolledApproximation : public readdy::model::actions::reactions::UncontrolledApproximation {
    using super = readdy::model::actions::reactions::UncontrolledApproximation;
public:
    CPUUncontrolledApproximation(CPUKernel* kernel, readdy::scalar timeStep);

    void perform(const util::PerformanceNode &node) override;

protected:
    CPUKernel *const kernel;
};
}
}
}
}
}
