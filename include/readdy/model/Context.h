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
 * Reactions and potentials come in two variants:
 *   - Internal, created by the responsible kernel
 *   - External, inserted from the 'outside', e.g. python prototypes
 *
 * The context registers both of them in separate maps. Before the simulation can start the
 * content of both of these maps is unified into a single map, which is referred to during the actual
 * run of the simulation.
 *
 * @file KernelContext.h
 * @brief Container class for time independent information of the KernelContext.
 * @author clonker
 * @author chrisfroe
 * @date 18.04.16
 * @todo write docs, is kbt really time indep (or should be treated as such)?
 */

#pragma once

#include <array>
#include <memory>
#include <vector>
#include <unordered_set>

#include <readdy/common/ParticleTypeTuple.h>
#include "ParticleTypeRegistry.h"

#include <readdy/model/potentials/PotentialRegistry.h>
#include <readdy/model/reactions/ReactionRegistry.h>
#include <readdy/model/topologies/TopologyRegistry.h>
#include <readdy/model/compartments/CompartmentRegistry.h>
#include <readdy/api/KernelConfiguration.h>

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

    const scalar &kBT() const;

    scalar &kBT();

    scalar boxVolume() const;

    const BoxSize &boxSize() const;

    BoxSize &boxSize();

    const PeriodicBoundaryConditions &periodicBoundaryConditions() const;

    PeriodicBoundaryConditions &periodicBoundaryConditions();

    std::tuple<Vec3, Vec3> getBoxBoundingVertices() const;

    const fix_pos_fun &fixPositionFun() const;

    const dist_squared_fun &distSquaredFun() const;

    const shortest_dist_fun &shortestDifferenceFun() const;

    const pbc_fun &applyPBCFun() const;

    const scalar calculateMaxCutoff() const;

    compartments::CompartmentRegistry &compartments();

    const compartments::CompartmentRegistry &compartments() const;

    /**
     * Copy the reactions and potentials of the internal and external registries into the actual registries, which
     * is used during the run of the simulation. This step is necessary before the simulation can start. Otherwise the
     * registered reactions and potentials will not take effect. The context is optionally logged in text format.
     * @param debugOutput decide if context information will be logged
     */
    void configure(bool debugOutput = false);

    /**
     * Returns whether reactions with positions shall be recorded in the state model, then obtainable by
     * the readdy::model::observables::Reactions observable.
     * @return whether reactions shall be recorded in the state model, by default false
     */
    const bool &recordReactionsWithPositions() const;

    /**
     * Returns whether reactions with positions shall be recorded in the state model, then obtainable by
     * the readdy::model::observables::Reactions observable.
     * @return whether reactions shall be recorded in the state model, by default false
     */
    bool &recordReactionsWithPositions();

    /**
     * Returns whether reaction counts shall be recorded in the state model (if it is supported). It is then obtainable
     * by the readdy::model::observables::ReactionCounts observable.
     * @return wheter reaction counts shall be recorded
     */
    const bool &recordReactionCounts() const;

    /**
     * Returns whether reaction counts shall be recorded in the state model (if it is supported). It is then obtainable
     * by the readdy::model::observables::ReactionCounts observable.
     * @return wheter reaction counts shall be recorded
     */
    bool &recordReactionCounts();

    top::TopologyRegistry &topology_registry();

    const top::TopologyRegistry &topology_registry() const;

    reactions::ReactionRegistry &reactions();

    const reactions::ReactionRegistry &reactions() const;

    ParticleTypeRegistry &particle_types();

    const ParticleTypeRegistry &particle_types() const;

    potentials::PotentialRegistry &potentials();

    const potentials::PotentialRegistry &potentials() const;

    KernelConfiguration &kernelConfiguration();

    const KernelConfiguration &kernelConfiguration() const;

    void setKernelConfiguration(const std::string &jsonStr);

    // ctor and dtor
    Context();

    ~Context();

    // move
    Context(Context &&rhs) = default;

    Context &operator=(Context &&rhs) = default;

    // copy
    Context(const Context &rhs) = default;

    Context &operator=(const Context &rhs) = default;

private:
    void updateFunctions();

    ParticleTypeRegistry _particleTypeRegistry;
    reactions::ReactionRegistry _reactionRegistry;
    potentials::PotentialRegistry _potentialRegistry;
    top::TopologyRegistry _topologyRegistry;
    compartments::CompartmentRegistry _compartmentRegistry;

    KernelConfiguration _kernelConfiguration;

    scalar _kBT{1};
    BoxSize _box_size{{1, 1, 1}};
    PeriodicBoundaryConditions _periodic_boundary{{true, true, true}};

    std::function<Vec3(const Vec3 &)> _pbc;
    std::function<void(Vec3 &)> _fixPositionFun;
    std::function<Vec3(const Vec3 &, const Vec3 &)> _diffFun;
    std::function<scalar(const Vec3 &, const Vec3 &)> _distFun;

    // here come horrible flags
    bool _recordReactionsWithPositions{false};
    bool _recordReactionCounts{false};

};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
