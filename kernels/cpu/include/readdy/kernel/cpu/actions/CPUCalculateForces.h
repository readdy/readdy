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
 * @date 14.07.16
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/common/thread/barrier.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
class CPUCalculateForces : public readdy::model::actions::CalculateForces {
    using data_bounds = std::tuple<data::EntryDataContainer::iterator, data::EntryDataContainer::iterator>;
    using nl_bounds = std::tuple<std::size_t, std::size_t>;
    using top_bounds = std::tuple<CPUStateModel::topologies_vec::const_iterator, CPUStateModel::topologies_vec::const_iterator>;
public:

    explicit CPUCalculateForces(CPUKernel *kernel) : kernel(kernel) {}

    void perform(const util::PerformanceNode &node) override;

protected:

    template<bool COMPUTE_VIRIAL>
    static void calculate_order2(std::size_t, nl_bounds nlBounds, CPUStateModel::data_type *data,
                                 const CPUStateModel::neighbor_list &nl, std::promise<scalar> &energyPromise,
                                 std::promise<Matrix33> &virialPromise,
                                 model::potentials::PotentialRegistry::potential_o2_registry pot2,
                                 model::Context::BoxSize box, model::Context::PeriodicBoundaryConditions pbc);

    static void calculate_topologies(std::size_t /*tid*/, top_bounds topBounds, model::top::TopologyActionFactory *taf,
                                     std::promise<scalar> &energyPromise);


    static void calculate_order1(std::size_t /*tid*/, data_bounds dataBounds,
                                 std::promise<scalar> &energyPromise, CPUStateModel::data_type *data,
                                 model::potentials::PotentialRegistry::potential_o1_registry pot1);

    CPUKernel *const kernel;
};
}
}
}
}
