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
 * This header file contains the definition of the ObservableFactory. Its purpose is to create observable of different
 * types in the form of unique_ptrs. The actual implementation of an observable can be changed by specializing the
 * dispatcher for its type and invoking a virtual (and then: overridden) method within the factory.
 *
 * @file ObservableFactory.h
 * @brief This header file contains the definition of the ObservableFactory.
 * @author clonker
 * @date 29.04.16
 */
#pragma once

#include <string>
#include <unordered_map>
#include <readdy/common/Utils.h>
#include <readdy/model/observables/Observables.h>
#include "Aggregators.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
class Kernel;

NAMESPACE_BEGIN(observables)

namespace detail {
template<typename T>
using is_observable_type = std::enable_if_t<std::is_base_of<model::observables::ObservableBase, T>::value>;
}

class ObservableFactory {
public:
    using stride_type = ObservableBase::stride_type;
    
    explicit ObservableFactory(Kernel *const kernel) : kernel(kernel) {};

    template<typename T>
    std::unique_ptr<Trivial<T>> collect(stride_type stride, T* observable, detail::is_observable_type<T>* = 0) const {
        return std::make_unique<Trivial<T>>(kernel, stride, observable);
    }
    
    std::unique_ptr<Energy> energy(stride_type stride) const {
        return std::make_unique<Energy>(kernel, stride);
    };

    virtual std::unique_ptr<Virial> virial(stride_type stride) const = 0;
    
    virtual std::unique_ptr<HistogramAlongAxis> histogramAlongAxis(stride_type stride, std::vector<scalar> binBorders, 
                                                                   std::vector<std::string> typesToCount, 
                                                                   unsigned int axis) const = 0;
    
    std::unique_ptr<NParticles> nParticles(stride_type stride) const { return nParticles(stride, {}); }
    
    virtual std::unique_ptr<NParticles> nParticles(stride_type stride, std::vector<std::string> typesToCount) const = 0;

    std::unique_ptr<Forces> forces(stride_type stride) const { return forces(stride, {}); }
    
    virtual std::unique_ptr<Forces> forces(stride_type stride, std::vector<std::string> typesToCount) const  = 0;

    std::unique_ptr<Positions> positions(stride_type stride) const { return positions(stride, {}); }

    virtual std::unique_ptr<Positions> positions(stride_type stride, std::vector<std::string> typesToCount) const = 0;

    virtual std::unique_ptr<RadialDistribution> radialDistribution(stride_type stride, std::vector<scalar> binBorders, 
                                                                   std::vector<std::string> typeCountFrom,
                                                                   std::vector<std::string> typeCountTo,
                                                                   scalar particleDensity) const = 0;

    virtual std::unique_ptr<Particles> particles(stride_type stride) const = 0;

    virtual std::unique_ptr<MeanSquaredDisplacement> msd(stride_type stride, std::vector<std::string> typesToCount, 
                                                         Particles *particlesObservable) const = 0;

    virtual std::unique_ptr<Reactions> reactions(stride_type stride) const = 0;

    virtual std::unique_ptr<ReactionCounts> reactionCounts(stride_type stride) const = 0;

    std::unique_ptr<Trajectory> trajectory(stride_type stride) const {
        return std::make_unique<Trajectory>(kernel, stride);
    }

    std::unique_ptr<FlatTrajectory> flatTrajectory(stride_type stride) const {
        return std::make_unique<FlatTrajectory>(kernel, stride);
    }

    std::unique_ptr<Topologies> topologies(stride_type stride) const {
        return std::make_unique<Topologies>(kernel, stride);
    }

protected:
    Kernel *const kernel;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
