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
 * @file SingleCPUObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */

#pragma once
#include <readdy/model/observables/ObservableFactory.h>

namespace readdy {
namespace kernel {
namespace scpu {
class SCPUKernel;
namespace observables {

class SCPUObservableFactory : public readdy::model::observables::ObservableFactory {

public:
    explicit SCPUObservableFactory(readdy::kernel::scpu::SCPUKernel* kernel);

    std::unique_ptr<readdy::model::observables::Virial> virial(stride_type stride) const override;

    std::unique_ptr<readdy::model::observables::HistogramAlongAxis>
    histogramAlongAxis(stride_type stride, std::vector<scalar> binBorders, std::vector<std::string> typesToCount,
                       unsigned int axis) const override;

    std::unique_ptr<readdy::model::observables::NParticles>
    nParticles(stride_type stride, std::vector<std::string> typesToCount) const override;

    std::unique_ptr<readdy::model::observables::Forces>
    forces(stride_type stride, std::vector<std::string> typesToCount) const override;

    std::unique_ptr<readdy::model::observables::Positions>
    positions(stride_type stride, std::vector<std::string> typesToCount) const override;

    std::unique_ptr<readdy::model::observables::RadialDistribution>
    radialDistribution(stride_type stride, std::vector<scalar> binBorders, std::vector<std::string> typeCountFrom,
                       std::vector<std::string> typeCountTo, scalar particleDensity) const override;

    std::unique_ptr<readdy::model::observables::Particles> particles(stride_type stride) const override;

    std::unique_ptr<readdy::model::observables::MeanSquaredDisplacement>
    msd(stride_type stride, std::vector<std::string> typesToCount,
        readdy::model::observables::Particles *particlesObservable) const override;

    std::unique_ptr<readdy::model::observables::Reactions> reactions(stride_type stride) const override;

    std::unique_ptr<readdy::model::observables::ReactionCounts> reactionCounts(stride_type stride) const override;
private:
    readdy::kernel::scpu::SCPUKernel *const kernel;
};

}
}
}
}
