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
 * The context is part of the simulation model. It contains all time-independent parameters, e.g. the size
 * and periodicity of the simulation box, definitions of particle species and the potentials
 * and reactions affecting them.
 *
 * Before the context object can be used in actions belonging to a kernel, the `configure` method has to be called
 * which will trigger certain rearrangements in the underlying registries.
 *
 * @file Context.h
 * @brief Container class for time independent information of the KernelContext.
 * @author clonker
 * @author chrisfroe
 * @date 18.04.16
 */

#pragma once

#include <array>
#include <memory>
#include <vector>
#include <unordered_set>

#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/api/KernelConfiguration.h>

#include "ParticleTypeRegistry.h"
#include <readdy/model/potentials/PotentialRegistry.h>
#include <readdy/model/reactions/ReactionRegistry.h>
#include <readdy/model/topologies/TopologyRegistry.h>
#include <readdy/model/compartments/CompartmentRegistry.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

class Context {
public:
    using BoxSize = std::array<scalar, 3>;
    using PeriodicBoundaryConditions = std::array<bool, 3>;
    using KernelConfiguration = conf::Configuration;

    using fix_pos_fun = std::function<void(Vec3 &)>;
    using pbc_fun = std::function<Vec3(const Vec3 &)>;
    using dist_squared_fun = std::function<scalar(const Vec3 &, const Vec3 &)>;
    using shortest_dist_fun = std::function<Vec3(const Vec3 &, const Vec3 &)>;

    const scalar &kBT() const {
        return _kBT;
    }

    scalar &kBT() {
        return _kBT;
    }

    scalar boxVolume() const {
        return _box_size.at(0) * _box_size.at(1) * _box_size.at(2);
    }

    const BoxSize &boxSize() const {
        return _box_size;
    }

    BoxSize &boxSize() {
        return _box_size;
    }

    const PeriodicBoundaryConditions &periodicBoundaryConditions() const {
        return _periodic_boundary;
    }

    PeriodicBoundaryConditions &periodicBoundaryConditions() {
        return _periodic_boundary;
    }

    std::tuple<Vec3, Vec3> getBoxBoundingVertices() const {
        const auto &boxSize = _box_size;
        Vec3 lowerLeft{static_cast<scalar>(-0.5) * boxSize[0],
                       static_cast<scalar>(-0.5) * boxSize[1],
                       static_cast<scalar>(-0.5) * boxSize[2]};
        auto upperRight = lowerLeft + Vec3(boxSize);
        return std::make_tuple(lowerLeft, upperRight);
    };

    const scalar calculateMaxCutoff() const;

    /**
     * Method for sanity checking some parts of the configuration.
     */
    void validate() const;

    /**
     * Construct a string that describes the context, i.e. particle types, reactions, potentials and topologies.
     */
    std::string describe();

    /**
     * Returns whether reactions with positions shall be recorded in the state model, then obtainable by
     * the readdy::model::observables::Reactions observable.
     * @return whether reactions shall be recorded in the state model, by default false
     */
    const bool &recordReactionsWithPositions() const {
        return _recordReactionsWithPositions;
    }

    /**
     * Returns whether reactions with positions shall be recorded in the state model, then obtainable by
     * the readdy::model::observables::Reactions observable.
     * @return whether reactions shall be recorded in the state model, by default false
     */
    bool &recordReactionsWithPositions() {
        return _recordReactionsWithPositions;
    }

    /**
     * Returns whether reaction counts shall be recorded in the state model (if it is supported). It is then obtainable
     * by the readdy::model::observables::ReactionCounts observable.
     * @return wheter reaction counts shall be recorded
     */
    const bool &recordReactionCounts() const {
        return _recordReactionCounts;
    }

    /**
     * Returns whether reaction counts shall be recorded in the state model (if it is supported). It is then obtainable
     * by the readdy::model::observables::ReactionCounts observable.
     * @return wheter reaction counts shall be recorded
     */
    bool &recordReactionCounts() {
        return _recordReactionCounts;
    }

    const bool &recordVirial() const {
        return _recordVirial;
    }

    bool &recordVirial() {
        return _recordVirial;
    }

    compartments::CompartmentRegistry &compartments() {
        return _compartmentRegistry;
    }

    const compartments::CompartmentRegistry &compartments() const {
        return _compartmentRegistry;
    }

    top::TopologyRegistry &topologyRegistry() {
        return _topologyRegistry;
    }

    const top::TopologyRegistry &topologyRegistry() const {
        return _topologyRegistry;
    }

    reactions::ReactionRegistry &reactions() {
        return _reactionRegistry;
    }

    const reactions::ReactionRegistry &reactions() const {
        return _reactionRegistry;
    }

    ParticleTypeRegistry &particleTypes() {
        return _particleTypeRegistry;
    }

    const ParticleTypeRegistry &particleTypes() const {
        return _particleTypeRegistry;
    }

    potentials::PotentialRegistry &potentials() {
        return _potentialRegistry;
    }

    const potentials::PotentialRegistry &potentials() const {
        return _potentialRegistry;
    }

    KernelConfiguration &kernelConfiguration() {
        return _kernelConfiguration;
    }

    const KernelConfiguration &kernelConfiguration() const {
        return _kernelConfiguration;
    }

    void setKernelConfiguration(const std::string &jsonStr);

    // ctor and dtor
    Context();

    ~Context() = default;

    // move
    Context(Context &&rhs) = default;

    Context &operator=(Context &&rhs) = default;

    // copy
    Context(const Context &rhs) = default;

    Context &operator=(const Context &rhs) = default;

private:
    ParticleTypeRegistry _particleTypeRegistry;
    reactions::ReactionRegistry _reactionRegistry;
    potentials::PotentialRegistry _potentialRegistry;
    top::TopologyRegistry _topologyRegistry;
    compartments::CompartmentRegistry _compartmentRegistry;

    KernelConfiguration _kernelConfiguration;

    scalar _kBT{1};
    BoxSize _box_size{{1, 1, 1}};
    PeriodicBoundaryConditions _periodic_boundary{{true, true, true}};

    // here come horrible flags
    bool _recordReactionsWithPositions{false};
    bool _recordReactionCounts{false};
    bool _recordVirial{false};

};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
