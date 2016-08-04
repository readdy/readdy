/**
 * todo
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
#include <readdy/model/_internal/ParticleTypePair.h>

#if BOOST_OS_MAC
#include <string>
#endif

namespace readdy {
    namespace model {
        class ParticleTypePairHasher {
        public:
            std::size_t operator()(const _internal::ParticleTypePair &k) const;
            std::size_t operator()(const std::tuple<unsigned int, unsigned int> &k) const;
        };

        class KernelContext {
        public:
            double getKBT() const;

            void setKBT(double kBT);

            std::array<double, 3> &getBoxSize() const;

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

            double getParticleRadius(const unsigned int &type) const;

            void setParticleRadius(const std::string &type, const double r);

            double getTimeStep() const;

            void setTimeStep(double dt);

            const std::vector<const reactions::Reaction<1> *> getAllOrder1Reactions() const;

            const reactions::Reaction<1> *const getReactionOrder1WithName(const std::string &name) const;

            const std::vector<reactions::Reaction<1>*> &getOrder1Reactions(const std::string &type) const;

            const std::vector<reactions::Reaction<1>*> &getOrder1Reactions(const unsigned int &type) const;

            const std::vector<const reactions::Reaction<2> *> getAllOrder2Reactions() const;

            const reactions::Reaction<2> *const getReactionOrder2WithName(const std::string &name) const;

            const std::vector<reactions::Reaction<2>*> &getOrder2Reactions(const std::string &type1, const std::string &type2) const;

            const std::vector<reactions::Reaction<2>*> &getOrder2Reactions(const unsigned int &type1, const unsigned int &type2) const;

            const boost::uuids::uuid &registerConversionReaction(const std::string &name, const std::string &from, const std::string &to, const double &rate);

            const boost::uuids::uuid &registerEnzymaticReaction(const std::string &name, const std::string &catalyst, const std::string &from, const std::string &to, const double &rate,
                                                                const double &eductDistance);

            const boost::uuids::uuid &registerFissionReaction(const std::string &name, const std::string &from, const std::string &to1, const std::string &to2, const double productDistance,
                                                              const double &rate, const double &weight1 = 0.5, const double &weight2 = 0.5);

            const boost::uuids::uuid &registerFusionReaction(const std::string &name, const std::string &from1, const std::string &from2, const std::string &to, const double &rate,
                                                             const double &eductDistance, const double &weight1 = 0.5, const double &weight2 = 0.5);

            const boost::uuids::uuid &registerDeathReaction(const std::string &name, const std::string &particleType, const double &rate);

            void deregisterPotential(const boost::uuids::uuid &potential);

            const boost::uuids::uuid &registerOrder1Potential(potentials::PotentialOrder1 const *const potential, const std::string &type);

            std::vector<potentials::PotentialOrder1*> getOrder1Potentials(const std::string &type) const;

            std::vector<potentials::PotentialOrder1*> getOrder1Potentials(const unsigned int type) const;

            const std::unordered_map<unsigned int, std::vector<potentials::PotentialOrder1*>> getAllOrder1Potentials() const;

            std::unordered_set<unsigned int> getAllOrder1RegisteredPotentialTypes() const;

            const boost::uuids::uuid &registerOrder2Potential(potentials::PotentialOrder2 const *const potential, const std::string &type1, const std::string &type2);

            const std::vector<potentials::PotentialOrder2*> &getOrder2Potentials(const std::string &type1, const std::string &type2) const;

            const std::vector<potentials::PotentialOrder2*> &getOrder2Potentials(const unsigned int type1, const unsigned int type2) const;

            const std::unordered_map<_internal::ParticleTypePair, std::vector<potentials::PotentialOrder2*>, readdy::model::ParticleTypePairHasher> getAllOrder2Potentials() const;

            std::vector<std::tuple<unsigned int, unsigned int>> getAllOrder2RegisteredPotentialTypes() const;

            std::vector<unsigned int> getAllRegisteredParticleTypes() const;
            std::string getParticleName(unsigned int id) const;

            void configure();

            // ctor and dtor
            KernelContext(reactions::ReactionFactory const *const reactionFactory);

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
        };

    }
}
#endif //READDY_MAIN_KERNELCONTEXT_H
