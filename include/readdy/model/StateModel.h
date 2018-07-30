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
 * The KernelStateModel keeps information about the current state of the system, like particle positions and forces.
 * A listener can be attached, that fires when the time step changes.
 *
 * @file KernelStateModel.h
 * @brief Defines the KernelStateModel, which gives information about the system's current state.
 * @author clonker
 * @date 18/04/16
 */

#pragma once
#include <vector>
#include <readdy/model/topologies/GraphTopology.h>
#include "Particle.h"
#include "readdy/common/ReaDDyVec3.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

class StateModel {
public:

    StateModel() = default;

    StateModel(const StateModel&) = delete;

    StateModel&operator=(const StateModel&) = delete;

    StateModel(StateModel&&) = default;

    StateModel& operator=(StateModel&&) = default;

    virtual ~StateModel() = default;

    // const accessor methods
    virtual const std::vector<Vec3> getParticlePositions() const = 0;

    virtual const std::vector<Particle> getParticles() const = 0;

    virtual Particle getParticleForIndex(std::size_t index) const = 0;

    virtual ParticleTypeId getParticleType(std::size_t index) const = 0;

    virtual void initializeNeighborList(scalar skin) = 0;

    virtual void updateNeighborList() = 0;

    virtual void clearNeighborList() = 0;

    virtual void addParticle(const Particle &p) = 0;

    virtual void addParticles(const std::vector<Particle> &p) = 0;

    virtual readdy::model::top::GraphTopology *const addTopology(TopologyTypeId type, const std::vector<TopologyParticle> &particles) = 0;

    virtual std::vector<Particle> getParticlesForTopology(const top::GraphTopology &topology) const;

    virtual std::vector<top::GraphTopology*> getTopologies() = 0;

    virtual top::GraphTopology const* getTopologyForParticle(top::Topology::particle_index particle) const = 0;

    virtual top::GraphTopology* getTopologyForParticle(top::Topology::particle_index particle) = 0;

    virtual void removeParticle(const Particle &p) = 0;

    virtual void removeAllParticles() = 0;

    virtual scalar energy() const = 0;

    virtual scalar &energy() = 0;

    virtual void toDenseParticleIndices(std::vector<std::size_t>::iterator begin,
                                        std::vector<std::size_t>::iterator end) const = 0;

    virtual void clear() = 0;
};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
