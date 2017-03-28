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
#include "Vec3.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

class KernelStateModel {
public:
    virtual ~KernelStateModel() = default;

    // const accessor methods
    virtual const std::vector<Vec3> getParticlePositions() const = 0;

    virtual const std::vector<Particle> getParticles() const = 0;

    virtual Particle getParticleForIndex(const std::size_t index) const = 0;

    virtual void updateNeighborList() = 0;

    virtual void clearNeighborList() = 0;

    virtual void calculateForces() = 0;

    virtual void addParticle(const Particle &p) = 0;

    virtual void addParticles(const std::vector<Particle> &p) = 0;

    virtual readdy::model::top::GraphTopology *const addTopology(const std::vector<TopologyParticle> &particles) = 0;

    virtual void removeParticle(const Particle &p) = 0;

    virtual void removeAllParticles() = 0;

    virtual double getEnergy() const = 0;

    virtual void expected_n_particles(const std::size_t n) = 0;
};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
