/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
                                 model::potentials::PotentialRegistry::PotentialsO2Map pot2,
                                 model::Context::BoxSize box, model::Context::PeriodicBoundaryConditions pbc);

    static void calculate_topologies(std::size_t /*tid*/, top_bounds topBounds, model::top::TopologyActionFactory *taf,
                                     std::promise<scalar> &energyPromise);


    static void calculate_order1(std::size_t /*tid*/, data_bounds dataBounds,
                                 std::promise<scalar> &energyPromise, CPUStateModel::data_type *data,
                                 model::potentials::PotentialRegistry::PotentialsO1Map pot1);

    CPUKernel *const kernel;
};
}
}
}
}
