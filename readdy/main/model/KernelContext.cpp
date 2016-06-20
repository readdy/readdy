/**
 * << detailed description >>
 *
 * @file KernelContext.cpp
 * @brief Implementation file of the KernelContext.
 * @author clonker
 * @date 18.04.16
 * @todo make proper reference to KernelContext.h, is kBT really indepdendent of t?
 */

#include <readdy/model/KernelContext.h>
#include <readdy/common/make_unique.h>
#include <unordered_map>
#include <boost/log/trivial.hpp>
#include <readdy/common/Utils.h>
#include <readdy/model/_internal/ParticleTypePair.h>
#include <readdy/model/reactions/Conversion.h>
#include <readdy/model/reactions/Enzymatic.h>
#include <readdy/model/reactions/Fission.h>
#include <readdy/model/reactions/Fusion.h>

namespace readdy {
    namespace model {
        struct ParticleTypePairHasher {
            std::size_t operator()(const _internal::ParticleTypePair &k) const {
                return hash_value(k);
            }

            std::size_t operator()(const std::tuple<unsigned int, unsigned int> &k) const {
                std::size_t seed = 0;
                const auto &t1 = std::get<0>(k);
                const auto &t2 = std::get<1>(k);
                if (t1 <= t2) {
                    boost::hash_combine(seed, t1);
                    boost::hash_combine(seed, t2);
                } else {
                    boost::hash_combine(seed, t2);
                    boost::hash_combine(seed, t1);
                }
                return seed;
            }
        };

        struct KernelContext::Impl {
            uint typeCounter;
            std::unordered_map<std::string, uint> typeMapping;
            double kBT = 0;
            std::array<double, 3> box_size{};
            std::array<bool, 3> periodic_boundary{};
            std::unordered_map<uint, double> diffusionConstants{};
            std::unordered_map<uint, double> particleRadii{};
            std::unordered_map<unsigned int, std::vector<std::unique_ptr<potentials::PotentialOrder1>>> potentialO1Registry{};
            std::unordered_map<_internal::ParticleTypePair, std::vector<std::unique_ptr<potentials::PotentialOrder2>>, ParticleTypePairHasher> potentialO2Registry{};

            std::unordered_map<unsigned int, std::vector<std::unique_ptr<reactions::Reaction>>> reactionOneEductRegistry{};
            std::unordered_map<_internal::ParticleTypePair, std::vector<std::unique_ptr<reactions::Reaction>>, ParticleTypePairHasher> reactionTwoEductsRegistry{};

            double timeStep;
        };


        double KernelContext::getKBT() const {
            return (*pimpl).kBT;
        }

        void KernelContext::setKBT(double kBT) {
            (*pimpl).kBT = kBT;
        }

        void KernelContext::setBoxSize(double dx, double dy, double dz) {
            (*pimpl).box_size = {dx, dy, dz};
        }

        void KernelContext::setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z) {
            (*pimpl).periodic_boundary = {pb_x, pb_y, pb_z};
        }

        KernelContext::KernelContext() : pimpl(std::make_unique<KernelContext::Impl>()) { }

        std::array<double, 3> &KernelContext::getBoxSize() const {
            return pimpl->box_size;
        }

        std::array<bool, 3> &KernelContext::getPeriodicBoundary() const {
            return pimpl->periodic_boundary;
        }

        double KernelContext::getDiffusionConstant(const std::string &particleType) const {
            return pimpl->diffusionConstants[pimpl->typeMapping[particleType]];
        }

        void KernelContext::setDiffusionConstant(const std::string &particleType, double D) {
            pimpl->diffusionConstants[getOrCreateTypeId(particleType)] = D;
        }

        double KernelContext::getTimeStep() const {
            return pimpl->timeStep;
        }

        void KernelContext::setTimeStep(double dt) {
            pimpl->timeStep = dt;
        }

        unsigned int KernelContext::getParticleTypeID(const std::string &name) const {
            return pimpl->typeMapping[name];
        }

        double KernelContext::getDiffusionConstant(uint particleType) const {
            return pimpl->diffusionConstants[particleType];
        }

        double KernelContext::getParticleRadius(const std::string &type) const {
            return getParticleRadius(pimpl->typeMapping[type]);
        }

        double KernelContext::getParticleRadius(const unsigned int &type) const {
            return pimpl->particleRadii[type];
        }

        void KernelContext::setParticleRadius(const std::string &particleType, const double r) {
            pimpl->particleRadii[getOrCreateTypeId(particleType)] = r;
        }

        unsigned int KernelContext::getOrCreateTypeId(const std::string &particleType) {
            uint t_id;
            if (pimpl->typeMapping.find(particleType) != pimpl->typeMapping.end()) {
                t_id = pimpl->typeMapping[particleType];
            } else {
                t_id = ++(pimpl->typeCounter);
                pimpl->typeMapping.emplace(particleType, t_id);
            }
            return t_id;
        }

        const boost::uuids::uuid &KernelContext::registerOrder2Potential(potentials::PotentialOrder2 const *const potential, const std::string &type1, const std::string &type2) {
            // wlog: type1 <= type2
            auto type1Id = pimpl->typeMapping[type1];
            auto type2Id = pimpl->typeMapping[type2];
            _internal::ParticleTypePair pp{type1Id, type2Id};
            if (pimpl->potentialO2Registry.find(pp) == pimpl->potentialO2Registry.end()) {
                pimpl->potentialO2Registry.emplace(pp, std::vector<std::unique_ptr<potentials::PotentialOrder2>>());
            }
            auto pot = potential->replicate();
            pot->configureForTypes(type1Id, type2Id);
            pimpl->potentialO2Registry[pp].push_back(std::unique_ptr<potentials::PotentialOrder2>(pot));
            return pimpl->potentialO2Registry[pp].back()->getId();
        }

        const std::vector<std::unique_ptr<potentials::PotentialOrder2>> &KernelContext::getOrder2Potentials(const std::string &type1, const std::string &type2) const {
            return getOrder2Potentials(pimpl->typeMapping[type1], pimpl->typeMapping[type2]);
        }

        const std::vector<std::unique_ptr<potentials::PotentialOrder2>> &KernelContext::getOrder2Potentials(const unsigned int type1, const unsigned int type2) const {
            return pimpl->potentialO2Registry[{type1, type2}];
        }

        std::unordered_set<std::tuple<unsigned int, unsigned int>, ParticleTypePairHasher> KernelContext::getAllOrder2RegisteredPotentialTypes() const {
            std::unordered_set<std::tuple<unsigned int, unsigned int>, ParticleTypePairHasher> result{};
            for (auto it = pimpl->potentialO2Registry.begin(); it != pimpl->potentialO2Registry.end(); ++it) {
                result.insert(std::make_tuple(it->first.t1, it->first.t2));
            }
            return result;
        }

        const boost::uuids::uuid &KernelContext::registerOrder1Potential(potentials::PotentialOrder1 const *const potential, const std::string &type) {
            auto typeId = pimpl->typeMapping[type];
            if (pimpl->potentialO1Registry.find(typeId) == pimpl->potentialO1Registry.end()) {
                pimpl->potentialO1Registry.insert(std::make_pair(typeId, std::vector<std::unique_ptr<potentials::PotentialOrder1>>()));
            }
            auto ptr = potential->replicate();
            ptr->configureForType(typeId);
            pimpl->potentialO1Registry[typeId].push_back(std::unique_ptr<potentials::PotentialOrder1>(ptr));
            return pimpl->potentialO1Registry[typeId].back()->getId();
        }

        const std::vector<std::unique_ptr<potentials::PotentialOrder1>> &KernelContext::getOrder1Potentials(const std::string &type) const {
            return getOrder1Potentials(pimpl->typeMapping[type]);
        }

        const std::vector<std::unique_ptr<potentials::PotentialOrder1>> &KernelContext::getOrder1Potentials(const unsigned int type) const {
            return pimpl->potentialO1Registry[type];
        }

        std::unordered_set<unsigned int> KernelContext::getAllOrder1RegisteredPotentialTypes() const {
            std::unordered_set<unsigned int> result{};
            for (auto it = pimpl->potentialO1Registry.begin(); it != pimpl->potentialO1Registry.end(); ++it) {
                result.insert(it->first);
            }
            return result;
        }

        void KernelContext::deregisterPotential(const boost::uuids::uuid &potential) {
            auto deleterO1 = [potential](std::unique_ptr<potentials::PotentialOrder1> const &p) -> bool { return potential == p->getId(); };
            auto deleterO2 = [potential](std::unique_ptr<potentials::PotentialOrder2> const &p) -> bool { return potential == p->getId(); };
            for (auto it = pimpl->potentialO1Registry.begin(); it != pimpl->potentialO1Registry.end(); ++it) {
                it->second.erase(std::remove_if(it->second.begin(), it->second.end(), deleterO1), it->second.end());
            }
            for (auto it = pimpl->potentialO2Registry.begin(); it != pimpl->potentialO2Registry.end(); ++it) {
                it->second.erase(std::remove_if(it->second.begin(), it->second.end(), deleterO2), it->second.end());
            }
        }

        const boost::uuids::uuid &KernelContext::registerConversionReaction(const std::string &name, const std::string &from, const std::string &to, const double &rate) {
            const auto &idFrom = pimpl->typeMapping[from];
            const auto &idTo = pimpl->typeMapping[to];
            if (pimpl->reactionOneEductRegistry.find(idFrom) == pimpl->reactionOneEductRegistry.end()) {
                pimpl->reactionOneEductRegistry.emplace(idFrom, std::vector<std::unique_ptr<reactions::Reaction>>());
            }
            pimpl->reactionOneEductRegistry[idFrom].push_back(std::make_unique<reactions::Conversion>(name, idFrom, idTo, rate));
            return pimpl->reactionOneEductRegistry[idFrom].back()->getId();
        }

        const boost::uuids::uuid &KernelContext::registerEnzymaticReaction(const std::string &name, const std::string &catalyst, const std::string &from, const std::string &to, const double &rate,
                                                                           const double &eductDistance) {
            const auto &idFrom = pimpl->typeMapping[from];
            const auto &idTo = pimpl->typeMapping[to];
            const auto &idCat = pimpl->typeMapping[catalyst];
            _internal::ParticleTypePair pp{idFrom, idCat};
            if (pimpl->reactionTwoEductsRegistry.find(pp) == pimpl->reactionTwoEductsRegistry.end()) {
                pimpl->reactionTwoEductsRegistry.emplace(pp, std::vector<std::unique_ptr<reactions::Reaction>>());
            }
            pimpl->reactionTwoEductsRegistry[pp].push_back(std::make_unique<reactions::Enzymatic>(name, idCat, idFrom, idTo, rate, eductDistance));
            return pimpl->reactionTwoEductsRegistry[pp].back()->getId();
        }

        const boost::uuids::uuid &KernelContext::registerFissionReaction(const std::string &name, const std::string &from, const std::string &to1, const std::string &to2, const double productDistance,
                                                                         const double &rate) {
            const auto &idFrom = pimpl->typeMapping[from];
            const auto &idTo1 = pimpl->typeMapping[to1];
            const auto &idTo2 = pimpl->typeMapping[to2];
            if (pimpl->reactionOneEductRegistry.find(idFrom) == pimpl->reactionOneEductRegistry.end()) {
                pimpl->reactionOneEductRegistry.emplace(idFrom, std::vector<std::unique_ptr<reactions::Reaction>>());
            }
            pimpl->reactionOneEductRegistry[idFrom].push_back(std::make_unique<reactions::Fission>(name, idFrom, idTo1, idTo2, productDistance, rate));
            return pimpl->reactionOneEductRegistry[idFrom].back()->getId();
        }

        const boost::uuids::uuid &KernelContext::registerFusionReaction(const std::string &name, const std::string &from1, const std::string &from2, const std::string &to, const double &rate,
                                                                        const double &eductDistance) {
            const auto &idFrom1 = pimpl->typeMapping[from1];
            const auto &idFrom2 = pimpl->typeMapping[from2];
            const auto &idTo = pimpl->typeMapping[to];
            _internal::ParticleTypePair pp{idFrom1, idFrom2};
            if (pimpl->reactionTwoEductsRegistry.find(pp) == pimpl->reactionTwoEductsRegistry.end()) {
                pimpl->reactionTwoEductsRegistry.emplace(pp, std::vector<std::unique_ptr<reactions::Reaction>>());
            }
            pimpl->reactionTwoEductsRegistry[pp].push_back(std::make_unique<reactions::Fusion>(name, idFrom1, idFrom2, idTo, rate, eductDistance));
            return pimpl->reactionTwoEductsRegistry[pp].back()->getId();
        }


        KernelContext &KernelContext::operator=(KernelContext &&rhs) = default;

        KernelContext::KernelContext(KernelContext &&rhs) = default;

        KernelContext::~KernelContext() = default;
    }
}






