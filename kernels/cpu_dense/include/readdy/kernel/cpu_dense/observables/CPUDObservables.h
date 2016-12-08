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
 * @date 22.11.16
 */

#ifndef READDY_KERNEL_CPU_DENSE_OBSERVABLES_H
#define READDY_KERNEL_CPU_DENSE_OBSERVABLES_H

#include <readdy/model/observables/Observables.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
class CPUDKernel;

namespace observables {

class CPUDPositions : public readdy::model::observables::Positions {
public:
    CPUDPositions(CPUDKernel *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {});

    virtual void evaluate() override;

protected:
    CPUDKernel *const kernel;
};

class CPUDParticles : public readdy::model::observables::Particles {
public:
    CPUDParticles(CPUDKernel *const kernel, unsigned int stride);

    virtual void evaluate() override;

protected:
    CPUDKernel *const kernel;
};

class CPUDHistogramAlongAxis : public readdy::model::observables::HistogramAlongAxis {

public:
    CPUDHistogramAlongAxis(CPUDKernel *const kernel, unsigned int stride,
                           const std::vector<double> &binBorders,
                           const std::vector<std::string> &typesToCount,
                           unsigned int axis);

    virtual void evaluate() override;

protected:
    CPUDKernel *const kernel;
    size_t size;
};

class CPUDNParticles : public readdy::model::observables::NParticles {
public:

    CPUDNParticles(CPUDKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});


    virtual void evaluate() override;

protected:
    CPUDKernel *const kernel;
};

class CPUDForces : public readdy::model::observables::Forces {
public:
    CPUDForces(CPUDKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {});

    virtual ~CPUDForces() {}

    virtual void evaluate() override;


protected:
    CPUDKernel *const kernel;
};


}
}
}
}
#endif //READDY_KERNEL_CPU_DENSE_OBSERVABLES_H
