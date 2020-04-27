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

namespace readdy::model {
class Kernel;
namespace observables {

namespace detail {
template<typename T>
using is_observable_type = std::enable_if_t<std::is_base_of<model::observables::ObservableBase, T>::value>;
}

class ObservableFactory {
protected:
    struct noop {
        template<typename T>
        void operator()(const T&) const {}
    };

    template<typename T>
    using ObsCallback = typename std::function<void(const typename T::result_type&)>;
public:

    explicit ObservableFactory(Kernel *const kernel) : kernel(kernel) {};

    [[nodiscard]] virtual std::unique_ptr<Energy> energy(Stride stride, ObsCallback<Energy> callback) const = 0;

    [[nodiscard]] std::unique_ptr<Energy> energy(Stride stride) const {
        return std::move(energy(stride, noop{}));
    }

    [[nodiscard]] virtual std::unique_ptr<Virial> virial(Stride stride, ObsCallback<Virial> callback) const = 0;

    [[nodiscard]] std::unique_ptr<Virial> virial(Stride stride) const {
        return std::move(virial(stride, noop{}));
    }

    [[nodiscard]] virtual std::unique_ptr<HistogramAlongAxis>
    histogramAlongAxis(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typesToCount,
                       unsigned int axis, ObsCallback<HistogramAlongAxis> callback) const = 0;

    [[nodiscard]] virtual std::unique_ptr<HistogramAlongAxis>
    histogramAlongAxis(Stride stride, const std::vector<scalar> &binBorders,
                       const std::vector<std::string> &typesToCount, unsigned int axis) const {
        return std::move(histogramAlongAxis(stride, binBorders, typesToCount, axis, noop{}));
    }
    
    [[nodiscard]] std::unique_ptr<NParticles>
    nParticles(Stride stride, ObsCallback<NParticles> callback) const {
        return std::move(nParticles(stride, {}, callback));
    }

    [[nodiscard]] std::unique_ptr<NParticles>
    nParticles(Stride stride, const std::vector<std::string> &typesToCount) const {
        return std::move(nParticles(stride, typesToCount, noop{}));
    }

    [[nodiscard]] std::unique_ptr<NParticles>
    nParticles(Stride stride) const {
        return std::move(nParticles(stride, {}, noop{}));
    }

    [[nodiscard]] virtual std::unique_ptr<NParticles>
    nParticles(Stride stride, std::vector<std::string> typesToCount, ObsCallback<NParticles> callback) const = 0;

    [[nodiscard]] std::unique_ptr<Forces>
    forces(Stride stride) const {
        return std::move(forces(stride, {}, noop{}));
    }

    [[nodiscard]] std::unique_ptr<Forces>
    forces(Stride stride, const std::vector<std::string> &typesToCount) const {
        return std::move(forces(stride, typesToCount, noop{}));
    }

    [[nodiscard]] std::unique_ptr<Forces>
    forces(Stride stride, const ObsCallback<Forces> &callback) const {
        return std::move(forces(stride, {}, callback));
    }

    [[nodiscard]] virtual std::unique_ptr<Forces>
    forces(Stride stride, std::vector<std::string> typesToCount, ObsCallback<Forces> callback) const  = 0;

    [[nodiscard]] std::unique_ptr<Positions>
    positions(Stride stride) const {
        return std::move(positions(stride, {}, noop{}));
    }

    [[nodiscard]] std::unique_ptr<Positions>
    positions(Stride stride, ObsCallback <Positions> callback) const {
        return std::move(positions(stride, {}, callback));
    }

    [[nodiscard]] std::unique_ptr<Positions>
    positions(Stride stride, std::vector<std::string> typesToCount) const {
        return std::move(positions(stride, typesToCount, noop{}));
    }

    [[nodiscard]] virtual std::unique_ptr<Positions>
    positions(Stride stride, std::vector<std::string> typesToCount, ObsCallback <Positions> callback) const = 0;

    [[nodiscard]] std::unique_ptr<RadialDistribution>
    radialDistribution(Stride stride, const std::vector<scalar> &binBorders, const std::vector<std::string> &typeCountFrom,
                       const std::vector<std::string> &typeCountTo, scalar particleDensity) const {
        return std::move(radialDistribution(stride, binBorders, typeCountFrom, typeCountTo, particleDensity, noop{}));
    }

    [[nodiscard]] virtual std::unique_ptr<RadialDistribution>
    radialDistribution(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typeCountFrom,
                       std::vector<std::string> typeCountTo, scalar particleDensity,
                       ObsCallback<RadialDistribution> callback) const = 0;

    [[nodiscard]] std::unique_ptr<Particles> particles(Stride stride) const {
        return std::move(particles(stride, noop{}));
    }

    [[nodiscard]] virtual std::unique_ptr<Particles> particles(Stride stride, ObsCallback<Particles> callback) const = 0;

    [[nodiscard]] std::unique_ptr<Reactions> reactions(Stride stride) const {
        return std::move(reactions(stride, noop{}));
    }

    [[nodiscard]] virtual std::unique_ptr<Reactions> reactions(Stride stride, ObsCallback<Reactions> callback) const = 0;

    [[nodiscard]] std::unique_ptr<ReactionCounts> reactionCounts(Stride stride) const {
        return std::move(reactionCounts(stride, noop{}));
    }

    [[nodiscard]] virtual std::unique_ptr<ReactionCounts> reactionCounts(Stride stride, ObsCallback<ReactionCounts> callback) const = 0;

    [[nodiscard]] std::unique_ptr<Trajectory> trajectory(Stride stride, ObsCallback<Trajectory> callback = [](const Trajectory::result_type&){}) const {
        auto obs = std::make_unique<Trajectory>(kernel, stride);
        obs->setCallback(callback);
        return std::move(obs);
    }

    [[nodiscard]] std::unique_ptr<FlatTrajectory> flatTrajectory(Stride stride, ObsCallback<FlatTrajectory> callback = [](const FlatTrajectory::result_type&){}) const {
        auto obs = std::make_unique<FlatTrajectory>(kernel, stride);
        obs->setCallback(callback);
        return std::move(obs);
    }

    [[nodiscard]] std::unique_ptr<Topologies> topologies(Stride stride, ObsCallback<Topologies> callback = [](const Topologies::result_type&){}) const {
        auto obs = std::make_unique<Topologies>(kernel, stride);
        obs->setCallback(callback);
        return std::move(obs);
    }

protected:
    Kernel *const kernel;
};

}
}
