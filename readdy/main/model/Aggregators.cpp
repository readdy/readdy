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
 * @file Aggregators.cpp
 * @brief Definition of several aggregators.
 * @author chrisfroe
 * @date 15.11.16
 */

#include <readdy/model/observables/Aggregators.h>

#include <readdy/model/Kernel.h>


namespace readdy {
namespace model {
namespace observables {

MeanSquaredDisplacement::MeanSquaredDisplacement(Kernel *const kernel, stride_type stride,
                                                 std::vector<std::string> typesToCount,
                                                 Particles *particlesObservable)
        : MeanSquaredDisplacement(kernel, stride, readdy::model::_internal::util::transformTypes2(typesToCount,
                                                                                                  kernel->context()),
                                  particlesObservable) {}

MeanSquaredDisplacement::MeanSquaredDisplacement(Kernel *const kernel, stride_type stride,
                                                 std::vector<ParticleTypeId> typesToCount,
                                                 Particles *particlesObservable)
        : Combiner(kernel, stride, particlesObservable), typesToCount(std::move(typesToCount)) {}

std::string MeanSquaredDisplacement::type() const {
    return "Mean squared displacement";
}

}
}
}
