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
 * @file CPUPotentialFactory.h
 * @brief This factory creates potentials for the cpu kernel.
 * @author clonker
 * @date 13.07.16
 */


#ifndef READDY_CPUKERNEL_CPUPOTENTIALFACTORY_H
#define READDY_CPUKERNEL_CPUPOTENTIALFACTORY_H

#include <readdy/model/potentials/PotentialFactory.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUKernel;
namespace potentials {
class CPUPotentialFactory : public readdy::model::potentials::PotentialFactory {
public:
    CPUPotentialFactory(CPUKernel *const kernel);
};
}
}
}
}
#endif //READDY_CPUKERNEL_CPUPOTENTIALFACTORY_H
