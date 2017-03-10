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
 * @file Observables.h
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#ifndef READDY_KERNEL_CPU_OBSERVABLES_H
#define READDY_KERNEL_CPU_OBSERVABLES_H

#include <readdy/model/observables/Observables.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUKernel;

namespace observables {

class CPUPositions : public readdy::model::observables::Positions {
public:
    CPUPositions(CPUKernel *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {});

    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUParticles : public readdy::model::observables::Particles {
public:
    CPUParticles(CPUKernel *const kernel, unsigned int stride);

    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUHistogramAlongAxis : public readdy::model::observables::HistogramAlongAxis {

public:
    CPUHistogramAlongAxis(CPUKernel *const kernel, unsigned int stride,
                       const std::vector<double> &binBorders,
                       const std::vector<std::string> &typesToCount,
                       unsigned int axis);

    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
    size_t size;
};

class CPUNParticles : public readdy::model::observables::NParticles {
public:

    CPUNParticles(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});


    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUForces : public readdy::model::observables::Forces {
public:
    CPUForces(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});

    virtual ~CPUForces() {}

    virtual void evaluate() override;


protected:
    CPUKernel *const kernel;
};

class CPUReactions : public readdy::model::observables::Reactions {
public:
    CPUReactions(CPUKernel *const kernel, unsigned int stride);

    virtual void evaluate() override;

protected:
    CPUKernel *const kernel;

};


}
}
}
}
#endif //READDY_KERNEL_CPU_OBSERVABLES_H
