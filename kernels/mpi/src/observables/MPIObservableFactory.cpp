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
 * @file MPIObservableFactory.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.05.19
 */

#include <readdy/kernel/mpi/observables/MPIObservableFactory.h>
#include <readdy/kernel/mpi/observables/MPIObservables.h>
#include <readdy/kernel/mpi/MPIKernel.h>

namespace readdy::kernel::mpi::observables {

MPIObservableFactory::MPIObservableFactory(MPIKernel *kernel) : readdy::model::observables::ObservableFactory(kernel),
                                                                kernel(kernel) {}

std::unique_ptr<readdy::model::observables::HistogramAlongAxis>
MPIObservableFactory::histogramAlongAxis(Stride stride, std::vector<scalar> binBorders,
                                         std::vector<std::string> typesToCount, unsigned int axis) const {
    return {std::make_unique<MPIHistogramAlongAxis>(kernel, stride, binBorders, typesToCount, axis)};
}

std::unique_ptr<readdy::model::observables::NParticles>
MPIObservableFactory::nParticles(Stride stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<MPINParticles>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::Forces>
MPIObservableFactory::forces(Stride stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<MPIForces>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::Positions>
MPIObservableFactory::positions(Stride stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<MPIPositions>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::RadialDistribution>
MPIObservableFactory::radialDistribution(Stride stride, std::vector<scalar> binBorders,
                                         std::vector<std::string> typeCountFrom, std::vector<std::string> typeCountTo,
                                         scalar particleDensity) const {
    return {std::make_unique<readdy::model::observables::RadialDistribution>(
            kernel, stride, binBorders, typeCountFrom, typeCountTo, particleDensity
    )};
}

std::unique_ptr<readdy::model::observables::Particles> MPIObservableFactory::particles(Stride stride) const {
    return {std::make_unique<MPIParticles>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::Reactions> MPIObservableFactory::reactions(Stride stride) const {
    return {std::make_unique<MPIReactions>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::ReactionCounts> MPIObservableFactory::reactionCounts(Stride stride) const {
    return {std::make_unique<MPIReactionCounts>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::Virial> MPIObservableFactory::virial(Stride stride) const {
    return {std::make_unique<MPIVirial>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::Energy> MPIObservableFactory::energy(Stride stride) const {
    return {std::make_unique<MPIEnergy>(kernel, stride)};
}

}
