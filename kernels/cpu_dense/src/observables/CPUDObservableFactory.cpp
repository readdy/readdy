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
 * @file ObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/observables/CPUDObservableFactory.h>
#include <readdy/kernel/cpu_dense/CPUDKernel.h>
#include <readdy/kernel/cpu_dense/observables/CPUDObservables.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace observables {
CPUDObservableFactory::CPUDObservableFactory(CPUDKernel *const kernel) : readdy::model::observables::ObservableFactory(kernel),
                                                             kernel(kernel) {
}

readdy::model::observables::NParticles *
CPUDObservableFactory::createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDNParticles(kernel, stride, typesToCount);
}

readdy::model::observables::HistogramAlongAxis *
CPUDObservableFactory::createAxisHistogramObservable(unsigned int stride, std::vector<double> binBorders,
                                                 std::vector<std::string> typesToCount,
                                                 unsigned int axis) const {
    return new CPUDHistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis);
}

readdy::model::observables::Forces *
CPUDObservableFactory::createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDForces(kernel, stride, typesToCount);
}

readdy::model::observables::ParticlePosition *
CPUDObservableFactory::createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDParticlePosition(kernel, stride, typesToCount);
}

}
}
}
}