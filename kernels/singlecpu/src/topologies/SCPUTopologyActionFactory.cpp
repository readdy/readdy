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
 * @file SCPUTopologyActionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/singlecpu/model/topologies/SCPUTopologyActionFactory.h>

namespace c_top = readdy::model::top;

using calc_harmonic_bond = c_top::CalculateHarmonicBondPotential;

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {
namespace top {

std::unique_ptr<calc_harmonic_bond> SCPUTopologyActionFactory::createCalculateHarmonicBondPotential() {
    return std::unique_ptr<calc_harmonic_bond>();
}
}
}
}
}
}
