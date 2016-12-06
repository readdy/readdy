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
 * @file ObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_DENSE_OBSERVABLEFACTORY_H
#define READDY_DENSE_OBSERVABLEFACTORY_H


#include <readdy/model/observables/ObservableFactory.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
class CPUDKernel;
namespace observables {
class CPUDObservableFactory : public readdy::model::observables::ObservableFactory {

public:
    CPUDObservableFactory(CPUDKernel *const kernel);

    virtual readdy::model::observables::NParticles *
    createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

    virtual readdy::model::observables::HistogramAlongAxis *
    createAxisHistogramObservable(unsigned int stride,
                                  std::vector<double> binBorders, std::vector<std::string> typesToCount,
                                  unsigned int axis) const override;

    virtual readdy::model::observables::Forces *
    createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

    virtual readdy::model::observables::ParticlePosition *
    createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

private:
    CPUDKernel *const kernel;
};
}
}
}
}

#endif //READDY_DENSE_OBSERVABLEFACTORY_H
