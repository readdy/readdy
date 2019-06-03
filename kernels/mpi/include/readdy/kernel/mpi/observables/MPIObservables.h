/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
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
 * « detailed description »
 *
 * @file MPIObservables.h
 * @brief « brief description »
 * @author chrisfroe
 * @date 03.06.19
 */

#pragma once


#include <readdy/model/observables/Observables.h>

namespace readdy::kernel::mpi {

class MPIKernel;

namespace observables {

class MPIVirial : public readdy::model::observables::Virial {
public:
    MPIVirial(MPIKernel *kernel, stride_type stride);

    void evaluate() override;

protected:
    MPIKernel * kernel;
};

class MPIPositions : public readdy::model::observables::Positions {
public:
    MPIPositions(MPIKernel* kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {});

    void evaluate() override;

protected:
    MPIKernel * kernel;
};

class MPIParticles : public readdy::model::observables::Particles {
public:
    MPIParticles(MPIKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    MPIKernel * kernel;
};

class MPIHistogramAlongAxis : public readdy::model::observables::HistogramAlongAxis {

public:
    MPIHistogramAlongAxis(MPIKernel* kernel, unsigned int stride,
                          const std::vector<scalar> &binBorders,
                          const std::vector<std::string> &typesToCount,
                          unsigned int axis);

    void evaluate() override;

protected:
    MPIKernel * kernel;
    size_t size;
};

class MPINParticles : public readdy::model::observables::NParticles {
public:

    MPINParticles(MPIKernel* kernel, unsigned int stride, std::vector<std::string> typesToCount = {});


    void evaluate() override;

protected:
    MPIKernel * kernel;
};

class MPIForces : public readdy::model::observables::Forces {
public:
    MPIForces(MPIKernel* kernel, unsigned int stride, std::vector<std::string> typesToCount = {});

    ~MPIForces() override = default;

    void evaluate() override;


protected:
    MPIKernel * kernel;
};

class MPIReactions : public readdy::model::observables::Reactions {
public:
    MPIReactions(MPIKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    MPIKernel * kernel;
};

class MPIReactionCounts : public readdy::model::observables::ReactionCounts {
public:
    MPIReactionCounts(MPIKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    MPIKernel * kernel;
};

}
}
