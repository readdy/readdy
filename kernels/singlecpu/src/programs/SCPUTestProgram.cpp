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
// Created by clonker on 11.04.16.
//

#include <readdy/kernel/singlecpu/programs/SCPUTestProgram.h>
#include <readdy/common/logging.h>

namespace sctp = readdy::kernel::scpu::programs;

sctp::SCPUTestProgram::SCPUTestProgram() : readdy::model::programs::Test() {

}

void sctp::SCPUTestProgram::execute() {
    log::console()->debug("execute called!");
}

sctp::SCPUTestProgram::~SCPUTestProgram() = default;







