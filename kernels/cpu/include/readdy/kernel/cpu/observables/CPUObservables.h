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
