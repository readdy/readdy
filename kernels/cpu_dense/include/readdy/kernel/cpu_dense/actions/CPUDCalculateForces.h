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
 * @file CalculateForces.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */


#ifndef READDY_DENSE_CALCULATEFORCES_H
#define READDY_DENSE_CALCULATEFORCES_H

#include <readdy/model/programs/Actions.h>
#include <readdy/kernel/cpu_dense/CPUDKernel.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace actions {
class CPUDCalculateForces : public readdy::model::actions::CalculateForces {

public:

    CPUDCalculateForces(CPUDKernel *const kernel) : kernel(kernel) {}

    virtual void perform() override {
        kernel->getKernelStateModel().calculateForces();
    }

protected:
    const CPUDKernel *const kernel;
};
}
}
}
}
#endif //READDY_DENSE_CALCULATEFORCES_H
