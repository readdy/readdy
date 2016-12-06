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
// Created by clonker on 08.04.16.
//

#ifndef READDY_MAIN_SINGLECPUPROGRAMFACTORY_H
#define READDY_MAIN_SINGLECPUPROGRAMFACTORY_H

#include <readdy/model/programs/ProgramFactory.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy {
namespace kernel {
namespace scpu {
class SCPUKernel;
namespace programs {
class SCPUProgramFactory : public readdy::model::programs::ProgramFactory {
public:
    SCPUProgramFactory(SCPUKernel *kernel);

private:
    SCPUKernel *kernel;
};
}
}
}
}


#endif //READDY_MAIN_SINGLECPUPROGRAMFACTORY_H
