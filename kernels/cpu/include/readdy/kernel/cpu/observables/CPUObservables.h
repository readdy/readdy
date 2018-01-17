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

#pragma once
#include <readdy/model/observables/Observables.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUKernel;

namespace observables {

class CPUVirial : public readdy::model::observables::Virial {
public:
    CPUVirial(CPUKernel *kernel, stride_type stride);

    void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUPositions : public readdy::model::observables::Positions {
public:
    CPUPositions(CPUKernel* kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {});

    void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUParticles : public readdy::model::observables::Particles {
public:
    CPUParticles(CPUKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUHistogramAlongAxis : public readdy::model::observables::HistogramAlongAxis {

public:
    CPUHistogramAlongAxis(CPUKernel* kernel, unsigned int stride,
                       const std::vector<scalar> &binBorders,
                       const std::vector<std::string> &typesToCount,
                       unsigned int axis);

    void evaluate() override;

protected:
    CPUKernel *const kernel;
    size_t size;
};

class CPUNParticles : public readdy::model::observables::NParticles {
public:

    CPUNParticles(CPUKernel* kernel, unsigned int stride, std::vector<std::string> typesToCount = {});


    void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUForces : public readdy::model::observables::Forces {
public:
    CPUForces(CPUKernel* kernel, unsigned int stride, std::vector<std::string> typesToCount = {});

    ~CPUForces() override = default;

    void evaluate() override;


protected:
    CPUKernel *const kernel;
};

class CPUReactions : public readdy::model::observables::Reactions {
public:
    CPUReactions(CPUKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    CPUKernel *const kernel;
};

class CPUReactionCounts : public readdy::model::observables::ReactionCounts {
public:
    CPUReactionCounts(CPUKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    CPUKernel *const kernel;
};

}
}
}
}
