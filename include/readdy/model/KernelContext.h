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

#ifndef READDY_MAIN_KERNELCONTEXT_H
#define READDY_MAIN_KERNELCONTEXT_H

#include <array>
#include <memory>
#include <vector>
#include <unordered_set>
#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/reactions/Reaction.h>
#include <readdy/model/reactions/ReactionFactory.h>
#include <readdy/common/ParticleTypePair.h>
#include <readdy/model/compartments/Compartment.h>

namespace readdy {
namespace model {

class UnknownParticleType : public std::runtime_error {
public:
    UnknownParticleType(const std::string &__arg) : runtime_error(__arg) {}
};

struct ParticleTypeInfo {
    std::string name;
    double diffusionConstant;
    double radius;
    readdy::model::Particle::flavor_t flavor;
    readdy::model::Particle::type_type typeId;

    ParticleTypeInfo(const std::string &name, const double diffusionConstant, const double radius,
                     const Particle::flavor_t flavor, const Particle::type_type typeId);
};

class KernelContext {
    using particle_t = readdy::model::Particle;

    using rea_ptr_vec1 = std::vector<std::unique_ptr<reactions::Reaction<1>>>;
    using rea_ptr_vec2 = std::vector<std::unique_ptr<reactions::Reaction<2>>>;

    using rdy_ptp = readdy::util::ParticleTypePair;
    using rdy_ptp_hasher = readdy::util::ParticleTypePairHasher;

    using reaction_o1_registry_internal = std::unordered_map<particle_t::type_type, rea_ptr_vec1>;
    using reaction_o2_registry_internal = std::unordered_map<rdy_ptp, rea_ptr_vec2, rdy_ptp_hasher>;

    using pot_ptr_vec1 = std::vector<std::unique_ptr<potentials::PotentialOrder1>>;
    using pot_ptr_vec1_external = std::vector<potentials::PotentialOrder1 *>;
    using pot_ptr_vec2 = std::vector<std::unique_ptr<potentials::PotentialOrder2>>;
    using pot_ptr_vec2_external = std::vector<potentials::PotentialOrder2 *>;

    using potential_o1_registry_internal = std::unordered_map<particle_t::type_type, pot_ptr_vec1>;
    using potential_o2_registry_internal = std::unordered_map<rdy_ptp, pot_ptr_vec2, rdy_ptp_hasher>;

public:

    using rdy_type_mapping = std::unordered_map<std::string, particle_t::type_type>;
    using rdy_reverse_type_mapping = std::unordered_map<particle_t::type_type, std::string>;

    using rdy_pot_1 = readdy::model::potentials::PotentialOrder1;
    using rdy_pot_1_registry = std::unordered_map<particle_t::type_type, std::vector<rdy_pot_1 *>>;

    using rdy_pot_2 = readdy::model::potentials::PotentialOrder2;
    using rdy_pot_2_registry = std::unordered_map<rdy_ptp, std::vector<rdy_pot_2 *>, rdy_ptp_hasher>;

    using reaction_o1_registry = std::unordered_map<particle_t::type_type, std::vector<reactions::Reaction<1> *>>;
    using reaction_o2_registry = std::unordered_map<rdy_ptp, std::vector<reactions::Reaction<2> *>, rdy_ptp_hasher>;

    using compartment_registry = std::vector<std::unique_ptr<readdy::model::compartments::Compartment>>;

    using fix_pos_fun = std::function<void(Vec3 &)>;
    using dist_squared_fun = std::function<double(const Vec3 &, const Vec3 &)>;
    using shortest_dist_fun = std::function<Vec3(const Vec3 &, const Vec3 &)>;

    double getKBT() const;

    void setKBT(double kBT);

    std::array<double, 3> &getBoxSize() const;

    std::tuple<readdy::model::Vec3, readdy::model::Vec3> getBoxBoundingVertices() const;

    void setBoxSize(double dx, double dy, double dz);

    const std::array<bool, 3> &getPeriodicBoundary() const;

    particle_t::type_type getParticleTypeID(const std::string &name) const;

    void setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z);

    const fix_pos_fun &getFixPositionFun() const;

    const dist_squared_fun &getDistSquaredFun() const;

    const shortest_dist_fun &getShortestDifferenceFun() const;

    void registerParticleType(const std::string &name, const double diffusionConst, const double radius,
                              const readdy::model::Particle::flavor_t flavor = readdy::model::Particle::FLAVOR_NORMAL);

    const ParticleTypeInfo& getParticleTypeInfo(const std::string& name) const;
    const ParticleTypeInfo& getParticleTypeInfo(const particle_t::type_type type) const;

    double getDiffusionConstant(const std::string &particleType) const;

    double getDiffusionConstant(const particle_t::type_type particleType) const;

    double getParticleRadius(const std::string &type) const;

    double getParticleRadius(const particle_t::type_type type) const;

    const std::vector<const reactions::Reaction<1> *> getAllOrder1Reactions() const;

    const reactions::Reaction<1> *const getReactionOrder1WithName(const std::string &name) const;

    const std::vector<reactions::Reaction<1> *> &getOrder1Reactions(const std::string &type) const;

    const std::vector<reactions::Reaction<1> *> &getOrder1Reactions(const particle_t::type_type type) const;

    const std::vector<const reactions::Reaction<2> *> getAllOrder2Reactions() const;

    const reactions::Reaction<2> *const getReactionOrder2WithName(const std::string &name) const;

    const std::vector<reactions::Reaction<2> *> &getOrder2Reactions(const std::string &type1,
                                                                    const std::string &type2) const;

    const std::vector<reactions::Reaction<2> *> &getOrder2Reactions(const particle_t::type_type type1,
                                                                    const particle_t::type_type type2) const;

    template<typename R>
    const short registerReaction(std::unique_ptr<R> r,
                                 typename std::enable_if<std::is_base_of<reactions::Reaction<1>, R>::value>::type * = 0) {
        log::console()->trace("registering reaction {}", *r);
        const auto id = r->getId();
        const auto type = r->getEducts()[0];
        if (reactionOneEductRegistryInternal.find(type) == reactionOneEductRegistryInternal.end()) {
            reactionOneEductRegistryInternal.emplace(type, rea_ptr_vec1());
        }
        reactionOneEductRegistryInternal[type].push_back(std::move(r));
        return id;
    }

    template<typename R>
    const short registerReaction(std::unique_ptr<R> r,
                                 typename std::enable_if<std::is_base_of<reactions::Reaction<2>, R>::value>::type * = 0) {
        log::console()->trace("registering reaction {}", *r);
        const auto id = r->getId();
        const auto t1 = r->getEducts()[0];
        const auto t2 = r->getEducts()[1];

        const readdy::util::ParticleTypePair pp{t1, t2};
        if (reactionTwoEductsRegistryInternal.find(pp) == reactionTwoEductsRegistryInternal.end()) {
            reactionTwoEductsRegistryInternal.emplace(pp, rea_ptr_vec2());
        }
        reactionTwoEductsRegistryInternal[pp].push_back(std::move(r));
        return id;
    }

    const short registerExternalReaction(reactions::Reaction<1> *r) {
        reactionOneEductRegistryExternal[r->getEducts()[0]].push_back(r);
        return r->getId();
    }

    const short registerExternalReaction(reactions::Reaction<2> *r) {
        reactionTwoEductsRegistryExternal[
                readdy::util::ParticleTypePair(r->getEducts()[0], r->getEducts()[1])
        ].push_back(r);
        return r->getId();
    }

    const short registerExternalPotential(potentials::PotentialOrder1 *potential) {
        auto typeId = getParticleTypeID(potential->particleType);
        if (potentialO1RegistryExternal.find(typeId) == potentialO1RegistryExternal.end()) {
            potentialO1RegistryExternal.emplace(std::make_pair(typeId, pot_ptr_vec1_external()));
        }
        potentialO1RegistryExternal[typeId].push_back(potential);
        return potential->getId();
    }

    const short registerExternalPotential(potentials::PotentialOrder2 *potential) {
        const auto id = potential->getId();
        auto type1Id = getParticleTypeID(potential->particleType1);
        auto type2Id = getParticleTypeID(potential->particleType2);
        readdy::util::ParticleTypePair pp{type1Id, type2Id};
        if (potentialO2RegistryExternal.find(pp) == potentialO2RegistryExternal.end()) {
            potentialO2RegistryExternal.emplace(pp, pot_ptr_vec2_external());
        }
        potentialO2RegistryExternal[pp].push_back(potential);
        return id;
    }

    template<typename R>
    const short registerPotential(std::unique_ptr<R> potential,
                                  typename std::enable_if<std::is_base_of<potentials::PotentialOrder1, R>::value>::type * = 0) {
        const auto id = potential->getId();
        auto typeId = getParticleTypeID(potential->particleType);
        if (potentialO1RegistryInternal.find(typeId) == potentialO1RegistryInternal.end()) {
            potentialO1RegistryInternal.insert(std::make_pair(typeId, pot_ptr_vec1()));
        }
        potentialO1RegistryInternal[typeId].push_back(std::move(potential));
        return id;
    }

    template<typename R>
    const short registerPotential(std::unique_ptr<R> potential,
                                  typename std::enable_if<std::is_base_of<potentials::PotentialOrder2, R>::value>::type * = 0) {
        const auto id = potential->getId();
        auto type1Id = getParticleTypeID(potential->particleType1);
        auto type2Id = getParticleTypeID(potential->particleType2);
        readdy::util::ParticleTypePair pp{type1Id, type2Id};
        if (potentialO2RegistryInternal.find(pp) == potentialO2RegistryInternal.end()) {
            potentialO2RegistryInternal.emplace(pp, pot_ptr_vec2());
        }
        potentialO2RegistryInternal[pp].push_back(std::move(potential));
        return id;
    }

    void deregisterPotential(const short);

    template<typename T>
    const short registerCompartment(std::unique_ptr<T> compartment) {
        // assert to prevent errors already at compile-time
        static_assert(std::is_base_of<compartments::Compartment, T>::value, "argument must be a compartment");
        const auto id = compartment->getId();
        compartmentRegistry->push_back(std::move(compartment));
        return id;
    };

    const compartment_registry &getCompartments() const;

    std::vector<rdy_pot_1 *> getOrder1Potentials(const std::string &type) const;

    std::vector<rdy_pot_1 *> getOrder1Potentials(const particle_t::type_type type) const;

    const rdy_pot_1_registry getAllOrder1Potentials() const;

    std::unordered_set<particle_t::type_type> getAllOrder1RegisteredPotentialTypes() const;

    const std::vector<rdy_pot_2 *> &getOrder2Potentials(const std::string &, const std::string &) const;

    const std::vector<rdy_pot_2 *> &getOrder2Potentials(const particle_t::type_type, const particle_t::type_type) const;

    const rdy_pot_2_registry getAllOrder2Potentials() const;

    std::vector<std::tuple<particle_t::type_type, particle_t::type_type>> getAllOrder2RegisteredPotentialTypes() const;

    /**
     * Get an unstructured list of all second order potentials. Useful in neighborlists, where maxcutoff is required.
     */
    const std::vector<potentials::PotentialOrder2 *> getVectorAllOrder2Potentials() const;

    std::vector<particle_t::type_type> getAllRegisteredParticleTypes() const;

    std::string getParticleName(particle_t::type_type id) const;

    /**
     * Copy the reactions and potentials of the internal and external registries into the actual registries, which
     * is used during the run of the simulation. This step is necessary before the simulation can start. Otherwise the
     * registered reactions and potentials will not take effect. The context is optionally logged in text format.
     * @param debugOutput decide if context information will be logged
     */
    void configure(bool debugOutput = false);

    const rdy_type_mapping &getTypeMapping() const;

    /**
     * Generate a map from particle_type_t to string. As there is no book-keeping of this reversed
     * structure, it is generated in place and then returned.
     */
    const rdy_reverse_type_mapping generateReverseTypeMapping() const;

    // ctor and dtor
    KernelContext();

    ~KernelContext();

    // move
    KernelContext(KernelContext &&rhs);

    KernelContext &operator=(KernelContext &&rhs);

    // copy
    KernelContext(const KernelContext &rhs) = delete;

    KernelContext &operator=(const KernelContext &rhs) = delete;

private:
    struct Impl;
    std::unique_ptr<readdy::model::KernelContext::Impl> pimpl;

    using reaction_o1_registry_external = reaction_o1_registry;
    using reaction_o2_registry_external = reaction_o2_registry;

    reaction_o1_registry reactionOneEductRegistry {};
    reaction_o1_registry_internal reactionOneEductRegistryInternal {};
    reaction_o1_registry_external reactionOneEductRegistryExternal {};
    reaction_o2_registry reactionTwoEductsRegistry {};
    reaction_o2_registry_internal reactionTwoEductsRegistryInternal {};
    reaction_o2_registry_external reactionTwoEductsRegistryExternal {};

    rdy_pot_1_registry potentialO1Registry {};
    rdy_pot_2_registry potentialO2Registry {};

    potential_o1_registry_internal potentialO1RegistryInternal {};
    rdy_pot_1_registry potentialO1RegistryExternal {};
    potential_o2_registry_internal potentialO2RegistryInternal {};
    rdy_pot_2_registry potentialO2RegistryExternal {};

    std::unique_ptr<compartment_registry> compartmentRegistry = std::make_unique<compartment_registry>();
};

}
}
#endif //READDY_MAIN_KERNELCONTEXT_H
