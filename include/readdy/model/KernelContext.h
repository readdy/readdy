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
#include <boost/predef.h>
#include <vector>
#include <unordered_set>
#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/reactions/Reaction.h>
#include <readdy/model/reactions/ReactionFactory.h>
#include <readdy/common/ParticleTypePair.h>

#if BOOST_OS_MAC
#include <string>
#endif

namespace readdy {
namespace model {

class KernelContext {
public:

    using rdy_type_mapping = std::unordered_map<std::string, uint>;

    using rdy_pot_1 = readdy::model::potentials::PotentialOrder1;
    using rdy_pot_1_registry = std::unordered_map<unsigned int, std::vector<rdy_pot_1 *>>;

    using rdy_pot_2 = readdy::model::potentials::PotentialOrder2;
    using rdy_pot_2_registry = std::unordered_map<readdy::util::ParticleTypePair,
            std::vector<rdy_pot_2 *>, readdy::util::ParticleTypePairHasher>;

    using reaction_o1_registry = std::unordered_map<unsigned int, std::vector<reactions::Reaction<1> *>>;
    using reaction_o2_registry = std::unordered_map<readdy::util::ParticleTypePair,
            std::vector<reactions::Reaction<2> *>, readdy::util::ParticleTypePairHasher>;

    double getKBT() const;

    void setKBT(double kBT);

    std::array<double, 3> &getBoxSize() const;

    std::tuple<readdy::model::Vec3, readdy::model::Vec3> getBoxBoundingVertices() const;

    void setBoxSize(double dx, double dy, double dz);

    const std::array<bool, 3> &getPeriodicBoundary() const;

    unsigned int getParticleTypeID(const std::string &name) const;

    void setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z);

    const std::function<void(Vec3 &)> &getFixPositionFun() const;

    const std::function<double(const Vec3 &, const Vec3 &)> &getDistSquaredFun() const;

    const std::function<Vec3(const Vec3 &, const Vec3 &)> &getShortestDifferenceFun() const;

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
    const boost::uuids::uuid &registerReaction(std::unique_ptr<R> r,
                                               typename std::enable_if<std::is_base_of<reactions::Reaction<1>, R>::value>::type * = 0) {
        reactionOneEductRegistryInternal->push_back(std::move(r));
        return registerReaction(r.get());
    }

    template<typename R>
    const boost::uuids::uuid &registerReaction(std::unique_ptr<R> r,
                                               typename std::enable_if<std::is_base_of<reactions::Reaction<2>, R>::value>::type * = 0) {
        reactionTwoEductsRegistryInternal->push_back(std::move(r));
        return registerReaction(r.get());
    }

    template<typename R>
    const boost::uuids::uuid &registerReaction(R* r,
                                               typename std::enable_if<std::is_base_of<reactions::Reaction<1>, R>::value>::type * = 0) {
        BOOST_LOG_TRIVIAL(trace) << "registering reaction " << *r;
        const auto type = r->getEducts()[0];
        if (reactionOneEductRegistry->find(type) == reactionOneEductRegistry->end()) {
            reactionOneEductRegistry->emplace(type, std::vector<reactions::Reaction<1>*>());
        }
        (*reactionOneEductRegistry)[type].push_back(r);
        return (*reactionOneEductRegistry)[type].back()->getId();
    };

    template<typename R>
    const boost::uuids::uuid &registerReaction(R* r,
                                               typename std::enable_if<std::is_base_of<reactions::Reaction<2>, R>::value>::type * = 0) {
        BOOST_LOG_TRIVIAL(trace) << "registering reaction " << *r;
        const auto t1 = r->getEducts()[0];
        const auto t2 = r->getEducts()[1];

        const readdy::util::ParticleTypePair pp{t1, t2};
        if (reactionTwoEductsRegistry->find(pp) == reactionTwoEductsRegistry->end()) {
            reactionTwoEductsRegistry->emplace(pp, std::vector<reactions::Reaction<2>*>());
        }
        (*reactionTwoEductsRegistry)[pp].push_back(r);
        return (*reactionTwoEductsRegistry)[pp].back()->getId();
    };

    template<typename R>
    const boost::uuids::uuid &registerPotential(R* potential, const std::string &type,
                                                typename std::enable_if<std::is_base_of<potentials::PotentialOrder1, R>::value>::type * = 0) {
        auto typeId = getOrCreateTypeId(type);
        if (potentialO1Registry->find(typeId) == potentialO1Registry->end()) {
            potentialO1Registry->insert(
                    std::make_pair(typeId, std::vector<potentials::PotentialOrder1*>()));
        }
        potential->configureForType(typeId);
        potentialO1Registry->operator[](typeId).push_back(std::move(potential));
        return potentialO1Registry->operator[](typeId).back()->getId();
    }

    template<typename R>
    const boost::uuids::uuid &registerPotential(R* potential,
                                                const std::string &type1, const std::string &type2,
                                                typename std::enable_if<std::is_base_of<potentials::PotentialOrder2, R>::value>::type * = 0) {
        // wlog: type1 <= type2
        auto type1Id = getOrCreateTypeId(type1);
        auto type2Id = getOrCreateTypeId(type2);
        readdy::util::ParticleTypePair pp{type1Id, type2Id};
        if (potentialO2Registry->find(pp) == potentialO2Registry->end()) {
            potentialO2Registry->emplace(pp, std::vector<potentials::PotentialOrder2*>());
        }
        potential->configureForTypes(type1Id, type2Id);
        potentialO2Registry->operator[](pp).push_back(std::move(potential));
        return potentialO2Registry->operator[](pp).back()->getId();
    }

    void deregisterPotential(const boost::uuids::uuid &potential);

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

    void configure();

    const rdy_type_mapping &getTypeMapping() const;

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
    using reaction_o1_registry_internal = std::vector<std::unique_ptr<reactions::Reaction<1>>>;
    using reaction_o2_registry_internal = std::vector<std::unique_ptr<reactions::Reaction<2>>>;

    using potential_o1_registry_internal = std::unordered_map<unsigned int, std::vector<std::unique_ptr<potentials::PotentialOrder1>>>;
    using potential_o2_registry_internal = std::unordered_map<readdy::util::ParticleTypePair, std::vector<std::unique_ptr<potentials::PotentialOrder2>>, readdy::util::ParticleTypePairHasher>;

    struct Impl;
    std::unique_ptr<readdy::model::KernelContext::Impl> pimpl;

    unsigned int getOrCreateTypeId(const std::string &name);

    std::unique_ptr<reaction_o1_registry> reactionOneEductRegistry = std::make_unique<reaction_o1_registry>();
    std::unique_ptr<reaction_o1_registry_internal> reactionOneEductRegistryInternal = std::make_unique<reaction_o1_registry_internal>();
    std::unique_ptr<reaction_o2_registry> reactionTwoEductsRegistry = std::make_unique<reaction_o2_registry>();
    std::unique_ptr<reaction_o2_registry_internal> reactionTwoEductsRegistryInternal = std::make_unique<reaction_o2_registry_internal>();

    std::unique_ptr<rdy_pot_1_registry> potentialO1Registry = std::make_unique<rdy_pot_1_registry>();
    std::unique_ptr<rdy_pot_2_registry> potentialO2Registry = std::make_unique<rdy_pot_2_registry>();
};

}
}
#endif //READDY_MAIN_KERNELCONTEXT_H
