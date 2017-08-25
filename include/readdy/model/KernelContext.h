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
#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/reactions/Reaction.h>
#include <readdy/model/reactions/ReactionFactory.h>
#include <readdy/model/reactions/ReactionRegistry.h>
#include <readdy/model/compartments/Compartment.h>
#include "Vec3.h"
#include "ParticleTypeRegistry.h"
#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/api/PotentialConfiguration.h>
#include <readdy/model/potentials/PotentialRegistry.h>
#include <readdy/model/topologies/TopologyTypeRegistry.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

class UnknownParticleType : public std::runtime_error {
public:
    explicit UnknownParticleType(const std::string &__arg) : runtime_error(__arg) {}
};

class KernelContext {
public:
    using compartment_registry = std::vector<std::unique_ptr<readdy::model::compartments::Compartment>>;

    using fix_pos_fun = std::function<void(Vec3 &)>;
    using pbc_fun = std::function<Vec3(const Vec3 &)>;
    using dist_squared_fun = std::function<scalar(const Vec3 &, const Vec3 &)>;
    using shortest_dist_fun = std::function<Vec3(const Vec3 &, const Vec3 &)>;

    scalar getKBT() const;

    void setKBT(scalar kBT);

    Vec3::data_arr &getBoxSize() const;

    std::tuple<readdy::model::Vec3, readdy::model::Vec3> getBoxBoundingVertices() const;

    void setBoxSize(scalar dx, scalar dy, scalar dz);

    const std::array<bool, 3> &getPeriodicBoundary() const;

    void setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z);

    const fix_pos_fun &getFixPositionFun() const;

    const dist_squared_fun &getDistSquaredFun() const;

    const shortest_dist_fun &getShortestDifferenceFun() const;

    const pbc_fun &getPBCFun() const;

    const scalar calculateMaxCutoff() const;

    template<typename T>
    const short registerCompartment(std::unique_ptr<T> compartment) {
        // assert to prevent errors already at compile-time
        static_assert(std::is_base_of<compartments::Compartment, T>::value, "argument must be a compartment");
        const auto id = compartment->getId();
        compartmentRegistry->push_back(std::move(compartment));
        return id;
    }

    const compartment_registry &getCompartments() const;

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

    void configureTopologyBondPotential(const std::string &type1, const std::string &type2, const api::Bond &bond);

    void configureTopologyAnglePotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                         const api::Angle &angle);

    void configureTopologyTorsionPotential(const std::string &type1, const std::string &type2, const std::string &type3,
                                           const std::string &type4, const api::TorsionAngle &torsionAngle);

    api::PotentialConfiguration &topology_potentials();

    const api::PotentialConfiguration &topology_potentials() const;

    top::TopologyTypeRegistry &topology_types();

    const top::TopologyTypeRegistry &topology_types() const;

    reactions::ReactionRegistry &reactions();

    const reactions::ReactionRegistry &reactions() const;

    ParticleTypeRegistry &particle_types();

    const ParticleTypeRegistry &particle_types() const;

    potentials::PotentialRegistry &potentials();

    const potentials::PotentialRegistry &potentials() const;

    // ctor and dtor
    KernelContext();

    ~KernelContext();

    // move
    KernelContext(KernelContext &&rhs) = delete;

    KernelContext &operator=(KernelContext &&rhs) = delete;

    // copy
    KernelContext(const KernelContext &rhs) = delete;

    KernelContext &operator=(const KernelContext &rhs) = delete;

private:
    struct Impl;
    std::unique_ptr<readdy::model::KernelContext::Impl> pimpl;

    api::PotentialConfiguration potentialConfiguration_{};

    ParticleTypeRegistry particleTypeRegistry_;
    reactions::ReactionRegistry reactionRegistry_;
    potentials::PotentialRegistry potentialRegistry_;
    top::TopologyTypeRegistry topologyTypes_;

    std::unique_ptr<compartment_registry> compartmentRegistry = std::make_unique<compartment_registry>();

    // here come horrible flags
    bool recordReactionsWithPositions_ = false;
    bool recordReactionCounts_ = false;
};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
