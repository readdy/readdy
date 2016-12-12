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
 * @file KernelContext.h
 * @brief Container class for time independent information of the KernelContext.
 * @author clonker
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

namespace readdy {
namespace model {

class UnknownParticleType : public std::runtime_error {
public:
    UnknownParticleType(const std::string &__arg) : runtime_error(__arg) {}
};

class KernelContext {
    using rea_ptr_vec1 = std::vector<std::unique_ptr<reactions::Reaction<1>>>;
    using rea_ptr_vec2 = std::vector<std::unique_ptr<reactions::Reaction<2>>>;

    using rdy_ptp = readdy::util::ParticleTypePair;
    using rdy_ptp_hasher = readdy::util::ParticleTypePairHasher;

    using reaction_o1_registry_internal = std::unordered_map<unsigned int, rea_ptr_vec1>;
    using reaction_o2_registry_internal = std::unordered_map<rdy_ptp, rea_ptr_vec2, rdy_ptp_hasher>;

    using pot_ptr_vec1 = std::vector<std::unique_ptr<potentials::PotentialOrder1>>;
    using pot_ptr_vec1_external = std::vector<potentials::PotentialOrder1 *>;
    using pot_ptr_vec2 = std::vector<std::unique_ptr<potentials::PotentialOrder2>>;
    using pot_ptr_vec2_external = std::vector<potentials::PotentialOrder2 *>;

    using potential_o1_registry_internal = std::unordered_map<unsigned int, pot_ptr_vec1>;
    using potential_o2_registry_internal = std::unordered_map<rdy_ptp, pot_ptr_vec2, rdy_ptp_hasher>;

public:

    using rdy_type_mapping = std::unordered_map<std::string, unsigned int>;
    using rdy_reverse_type_mapping = std::unordered_map<unsigned int, std::string>;

    using rdy_pot_1 = readdy::model::potentials::PotentialOrder1;
    using rdy_pot_1_registry = std::unordered_map<unsigned int, std::vector<rdy_pot_1 *>>;

    using rdy_pot_2 = readdy::model::potentials::PotentialOrder2;
    using rdy_pot_2_registry = std::unordered_map<rdy_ptp, std::vector<rdy_pot_2 *>, rdy_ptp_hasher>;

    using reaction_o1_registry = std::unordered_map<unsigned int, std::vector<reactions::Reaction<1> *>>;
    using reaction_o2_registry = std::unordered_map<rdy_ptp, std::vector<reactions::Reaction<2> *>, rdy_ptp_hasher>;

    using fix_pos_fun = std::function<void(Vec3 &)>;
    using dist_squared_fun = std::function<double(const Vec3 &, const Vec3 &)>;
    using shortest_dist_fun = std::function<Vec3(const Vec3 &, const Vec3 &)>;

    double getKBT() const;

    void setKBT(double kBT);

    std::array<double, 3> &getBoxSize() const;

    std::tuple<readdy::model::Vec3, readdy::model::Vec3> getBoxBoundingVertices() const;

    void setBoxSize(double dx, double dy, double dz);

    const std::array<bool, 3> &getPeriodicBoundary() const;

    unsigned int getParticleTypeID(const std::string &name) const;

    void setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z);

    const fix_pos_fun &getFixPositionFun() const;

    const dist_squared_fun &getDistSquaredFun() const;

    const shortest_dist_fun &getShortestDifferenceFun() const;

    double getDiffusionConstant(const std::string &particleType) const;

    double getDiffusionConstant(const unsigned int particleType) const;

    void setDiffusionConstant(const std::string &particleType, double D);

    double getParticleRadius(const std::string &type) const;

    double getParticleRadius(const unsigned int type) const;

    void setParticleRadius(const std::string &type, const double r);

    double getTimeStep() const;

    void setTimeStep(double dt);

    const std::vector<const reactions::Reaction<1> *> getAllOrder1Reactions() const;

    const reactions::Reaction<1> *const getReactionOrder1WithName(const std::string &name) const;

    const std::vector<reactions::Reaction<1> *> &getOrder1Reactions(const std::string &type) const;

    const std::vector<reactions::Reaction<1> *> &getOrder1Reactions(const unsigned int type) const;

    const std::vector<const reactions::Reaction<2> *> getAllOrder2Reactions() const;

    const reactions::Reaction<2> *const getReactionOrder2WithName(const std::string &name) const;

    const std::vector<reactions::Reaction<2> *> &getOrder2Reactions(const std::string &type1,
                                                                    const std::string &type2) const;

    const std::vector<reactions::Reaction<2> *> &getOrder2Reactions(const unsigned int type1,
                                                                    const unsigned int type2) const;

    template<typename R>
    const short registerReaction(std::unique_ptr<R> r,
                                 typename std::enable_if<std::is_base_of<reactions::Reaction<1>, R>::value>::type * = 0) {
        log::console()->trace("registering reaction {}", *r);
        const auto id = r->getId();
        const auto type = r->getEducts()[0];
        if (reactionOneEductRegistryInternal->find(type) == reactionOneEductRegistryInternal->end()) {
            reactionOneEductRegistryInternal->emplace(type, rea_ptr_vec1());
        }
        (*reactionOneEductRegistryInternal)[type].push_back(std::move(r));
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
        if (reactionTwoEductsRegistryInternal->find(pp) == reactionTwoEductsRegistryInternal->end()) {
            reactionTwoEductsRegistryInternal->emplace(pp, rea_ptr_vec2());
        }
        (*reactionTwoEductsRegistryInternal)[pp].push_back(std::move(r));
        return id;
    }

    const short registerExternalReaction(reactions::Reaction<1> *r) {
        (*reactionOneEductRegistryExternal)[r->getEducts()[0]].push_back(r);
        return r->getId();
    }

    const short registerExternalReaction(reactions::Reaction<2> *r) {
        (*reactionTwoEductsRegistryExternal)[
                readdy::util::ParticleTypePair(r->getEducts()[0], r->getEducts()[1])
        ].push_back(r);
        return r->getId();
    }

    const short registerExternalPotential(potentials::PotentialOrder1 *potential, const std::string &type) {
        auto typeId = getOrCreateTypeId(type);
        if (potentialO1RegistryExternal->find(typeId) == potentialO1RegistryExternal->end()) {
            potentialO1RegistryExternal->emplace(std::make_pair(typeId, pot_ptr_vec1_external()));
        }
        (*potentialO1RegistryExternal)[typeId].push_back(potential);
        return potential->getId();
    }

    const short registerExternalPotential(potentials::PotentialOrder2 *potential,
                                          const std::string &type1, const std::string &type2) {
        const auto id = potential->getId();
        auto type1Id = getOrCreateTypeId(type1);
        auto type2Id = getOrCreateTypeId(type2);
        readdy::util::ParticleTypePair pp{type1Id, type2Id};
        if (potentialO2RegistryExternal->find(pp) == potentialO2RegistryExternal->end()) {
            potentialO2RegistryExternal->emplace(pp, pot_ptr_vec2_external());
        }
        (*potentialO2RegistryExternal)[pp].push_back(potential);
        return id;
    }

    template<typename R>
    const short registerPotential(std::unique_ptr<R> potential, const std::string &type,
                                  typename std::enable_if<std::is_base_of<potentials::PotentialOrder1, R>::value>::type * = 0) {
        const auto id = potential->getId();
        auto typeId = getOrCreateTypeId(type);
        if (potentialO1RegistryInternal->find(typeId) == potentialO1RegistryInternal->end()) {
            potentialO1RegistryInternal->insert(std::make_pair(typeId, pot_ptr_vec1()));
        }
        (*potentialO1RegistryInternal)[typeId].push_back(std::move(potential));
        return id;
    }

    template<typename R>
    const short registerPotential(std::unique_ptr<R> potential, const std::string &type1, const std::string &type2,
                                  typename std::enable_if<std::is_base_of<potentials::PotentialOrder2, R>::value>::type * = 0) {
        const auto id = potential->getId();
        auto type1Id = getOrCreateTypeId(type1);
        auto type2Id = getOrCreateTypeId(type2);
        readdy::util::ParticleTypePair pp{type1Id, type2Id};
        if (potentialO2RegistryInternal->find(pp) == potentialO2RegistryInternal->end()) {
            potentialO2RegistryInternal->emplace(pp, pot_ptr_vec2());
        }
        (*potentialO2RegistryInternal)[pp].push_back(std::move(potential));
        return id;
    }

    void deregisterPotential(const short);

    std::vector<rdy_pot_1 *> getOrder1Potentials(const std::string &type) const;

    std::vector<rdy_pot_1 *> getOrder1Potentials(const unsigned int type) const;

    const rdy_pot_1_registry getAllOrder1Potentials() const;

    std::unordered_set<unsigned int> getAllOrder1RegisteredPotentialTypes() const;

    const std::vector<rdy_pot_2 *> &getOrder2Potentials(const std::string &, const std::string &) const;

    const std::vector<rdy_pot_2 *> &getOrder2Potentials(const unsigned int, const unsigned int) const;

    const rdy_pot_2_registry getAllOrder2Potentials() const;

    std::vector<std::tuple<unsigned int, unsigned int>> getAllOrder2RegisteredPotentialTypes() const;

    /**
     * Get an unstructured list of all second order potentials. Useful in neighborlists, where maxcutoff is required.
     */
    const std::vector<potentials::PotentialOrder2 *> getVectorAllOrder2Potentials() const;

    std::vector<unsigned int> getAllRegisteredParticleTypes() const;

    std::string getParticleName(unsigned int id) const;

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

    unsigned int getOrCreateTypeId(const std::string &name);

    using reaction_o1_registry_external = reaction_o1_registry;
    using reaction_o2_registry_external = reaction_o2_registry;

    std::unique_ptr<reaction_o1_registry> reactionOneEductRegistry = std::make_unique<reaction_o1_registry>();
    std::unique_ptr<reaction_o1_registry_internal> reactionOneEductRegistryInternal = std::make_unique<reaction_o1_registry_internal>();
    std::unique_ptr<reaction_o1_registry_external> reactionOneEductRegistryExternal = std::make_unique<reaction_o1_registry_external>();
    std::unique_ptr<reaction_o2_registry> reactionTwoEductsRegistry = std::make_unique<reaction_o2_registry>();
    std::unique_ptr<reaction_o2_registry_internal> reactionTwoEductsRegistryInternal = std::make_unique<reaction_o2_registry_internal>();
    std::unique_ptr<reaction_o2_registry_external> reactionTwoEductsRegistryExternal = std::make_unique<reaction_o2_registry_external>();

    std::unique_ptr<rdy_pot_1_registry> potentialO1Registry = std::make_unique<rdy_pot_1_registry>();
    std::unique_ptr<rdy_pot_2_registry> potentialO2Registry = std::make_unique<rdy_pot_2_registry>();

    std::unique_ptr<potential_o1_registry_internal> potentialO1RegistryInternal = std::make_unique<potential_o1_registry_internal>();
    std::unique_ptr<rdy_pot_1_registry> potentialO1RegistryExternal = std::make_unique<rdy_pot_1_registry>();
    std::unique_ptr<potential_o2_registry_internal> potentialO2RegistryInternal = std::make_unique<potential_o2_registry_internal>();
    std::unique_ptr<rdy_pot_2_registry> potentialO2RegistryExternal = std::make_unique<rdy_pot_2_registry>();
};

}
}
#endif //READDY_MAIN_KERNELCONTEXT_H
