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
#include <readdy/common/Utils.h>

namespace readdy {
    namespace model {
        struct KernelContext::Impl {
            uint typeCounter;
            std::unordered_map<std::string, unsigned int> typeMapping;
            double kBT = 1;
            std::array<double, 3> box_size{{1, 1, 1}};
            std::array<bool, 3> periodic_boundary{{true, true, true}};
            std::unordered_map<uint, double> diffusionConstants{};
            std::unordered_map<uint, double> particleRadii{};

            std::unordered_map<unsigned int, std::vector<potentials::PotentialOrder1 *>> potentialO1Registry{};
            std::unordered_map<readdy::util::ParticleTypePair, std::vector<potentials::PotentialOrder2 *>, readdy::util::ParticleTypePairHasher> potentialO2Registry{};
            std::unordered_map<unsigned int, std::vector<std::unique_ptr<potentials::PotentialOrder1>>> internalPotentialO1Registry{};
            std::unordered_map<readdy::util::ParticleTypePair, std::vector<std::unique_ptr<potentials::PotentialOrder2>>, readdy::util::ParticleTypePairHasher> internalPotentialO2Registry{};

            double timeStep;

            std::function<void(Vec3 &)> fixPositionFun = [](
                    Vec3 &vec) -> void { readdy::model::fixPosition<true, true, true>(vec, 1., 1., 1.); };
            std::function<Vec3(const Vec3 &, const Vec3 &)> diffFun = [](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                return readdy::model::shortestDifference<true, true, true>(lhs, rhs, 1., 1., 1.);
            };
            std::function<double(const Vec3 &, const Vec3 &)> distFun = [&](const Vec3 &lhs,
                                                                            const Vec3 &rhs) -> double {
                auto dv = diffFun(lhs, rhs);
                return dv * dv;
            };

            std::vector<reactions::Reaction<1> *> defaultReactionsO1 {};
            std::vector<reactions::Reaction<2> *> defaultReactionsO2 {};
            std::vector<potentials::PotentialOrder1 *> defaultPotentialsO1 {};
            std::vector<potentials::PotentialOrder2 *> defaultPotentialsO2 {};

            void updateDistAndFixPositionFun() {
                if (periodic_boundary[0]) {
                    if (periodic_boundary[1]) {
                        if (periodic_boundary[2]) {
                            diffFun = [&](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                                return readdy::model::shortestDifference<true, true, true>(lhs, rhs, box_size[0],
                                                                                           box_size[1], box_size[2]);
                            };
                            fixPositionFun = [&](Vec3 &vec) -> void {
                                readdy::model::fixPosition<true, true, true>(vec, box_size[0], box_size[1],
                                                                             box_size[2]);
                            };
                        } else {
                            diffFun = [&](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                                return readdy::model::shortestDifference<true, true, false>(lhs, rhs, box_size[0],
                                                                                            box_size[1], box_size[2]);
                            };
                            fixPositionFun = [&](Vec3 &vec) -> void {
                                readdy::model::fixPosition<true, true, false>(vec, box_size[0], box_size[1],
                                                                              box_size[2]);
                            };
                        }
                    } else {
                        if (periodic_boundary[2]) {
                            diffFun = [&](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                                return readdy::model::shortestDifference<true, false, true>(lhs, rhs, box_size[0],
                                                                                            box_size[1], box_size[2]);
                            };
                            fixPositionFun = [&](Vec3 &vec) -> void {
                                readdy::model::fixPosition<true, false, true>(vec, box_size[0], box_size[1],
                                                                              box_size[2]);
                            };
                        } else {
                            diffFun = [&](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                                return readdy::model::shortestDifference<true, false, false>(lhs, rhs, box_size[0],
                                                                                             box_size[1], box_size[2]);
                            };
                            fixPositionFun = [&](Vec3 &vec) -> void {
                                readdy::model::fixPosition<true, false, false>(vec, box_size[0], box_size[1],
                                                                               box_size[2]);
                            };
                        }
                    }
                } else {
                    if (periodic_boundary[1]) {
                        if (periodic_boundary[2]) {
                            diffFun = [&](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                                return readdy::model::shortestDifference<false, true, true>(lhs, rhs, box_size[0],
                                                                                            box_size[1], box_size[2]);
                            };
                            fixPositionFun = [&](Vec3 &vec) -> void {
                                readdy::model::fixPosition<false, true, true>(vec, box_size[0], box_size[1],
                                                                              box_size[2]);
                            };
                        } else {
                            diffFun = [&](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                                return readdy::model::shortestDifference<false, true, false>(lhs, rhs, box_size[0],
                                                                                             box_size[1], box_size[2]);
                            };
                            fixPositionFun = [&](Vec3 &vec) -> void {
                                readdy::model::fixPosition<false, true, false>(vec, box_size[0], box_size[1],
                                                                               box_size[2]);
                            };
                        }
                    } else {
                        if (periodic_boundary[2]) {
                            diffFun = [&](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                                return readdy::model::shortestDifference<false, false, true>(lhs, rhs, box_size[0],
                                                                                             box_size[1], box_size[2]);
                            };
                            fixPositionFun = [&](Vec3 &vec) -> void {
                                readdy::model::fixPosition<false, false, true>(vec, box_size[0], box_size[1],
                                                                               box_size[2]);
                            };
                        } else {
                            diffFun = [&](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
                                return readdy::model::shortestDifference<false, false, false>(lhs, rhs, box_size[0],
                                                                                              box_size[1], box_size[2]);
                            };
                            fixPositionFun = [&](Vec3 &vec) -> void {
                                readdy::model::fixPosition<false, false, false>(vec, box_size[0], box_size[1],
                                                                                box_size[2]);
                            };
                        }
                    }
                }
            }
        };


        double KernelContext::getKBT() const {
            return (*pimpl).kBT;
        }

        void KernelContext::setKBT(double kBT) {
            (*pimpl).kBT = kBT;
        }

        void KernelContext::setBoxSize(double dx, double dy, double dz) {
            (*pimpl).box_size = {dx, dy, dz};
            pimpl->updateDistAndFixPositionFun();
        }

        void KernelContext::setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z) {
            (*pimpl).periodic_boundary = {pb_x, pb_y, pb_z};
            pimpl->updateDistAndFixPositionFun();
        }

        KernelContext::KernelContext() : pimpl(std::make_unique<KernelContext::Impl>()) {
        }

        std::array<double, 3> &KernelContext::getBoxSize() const {
            return pimpl->box_size;
        }

        const std::array<bool, 3> &KernelContext::getPeriodicBoundary() const {
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

        double KernelContext::getParticleRadius(const unsigned int type) const {
            if (pimpl->particleRadii.find(type) == pimpl->particleRadii.end()) {
                BOOST_LOG_TRIVIAL(warning) << "No particle radius was set for the particle type id " << type
                                           << ", setting r=1";
                pimpl->particleRadii[type] = 1;
            }
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

        const boost::uuids::uuid &
        KernelContext::registerOrder2Potential(potentials::PotentialOrder2 const *const potential,
                                               const std::string &type1, const std::string &type2) {
            // wlog: type1 <= type2
            auto type1Id = pimpl->typeMapping[type1];
            auto type2Id = pimpl->typeMapping[type2];
            readdy::util::ParticleTypePair pp{type1Id, type2Id};
            if (pimpl->internalPotentialO2Registry.find(pp) == pimpl->internalPotentialO2Registry.end()) {
                pimpl->internalPotentialO2Registry.emplace(pp,
                                                           std::vector<std::unique_ptr<potentials::PotentialOrder2>>());
            }
            auto pot = potential->replicate();
            pot->configureForTypes(type1Id, type2Id);
            pimpl->internalPotentialO2Registry[pp].push_back(std::unique_ptr<potentials::PotentialOrder2>(pot));
            return pimpl->internalPotentialO2Registry[pp].back()->getId();
        }

        const std::vector<potentials::PotentialOrder2 *> &
        KernelContext::getOrder2Potentials(const std::string &type1, const std::string &type2) const {
            return getOrder2Potentials(pimpl->typeMapping[type1], pimpl->typeMapping[type2]);
        }

        const std::vector<potentials::PotentialOrder2 *> &
        KernelContext::getOrder2Potentials(const unsigned int type1, const unsigned int type2) const {
            return readdy::utils::collections::getOrDefault(pimpl->potentialO2Registry, {type1, type2}, pimpl->defaultPotentialsO2);
        }

        std::vector<std::tuple<unsigned int, unsigned int>>
        KernelContext::getAllOrder2RegisteredPotentialTypes() const {
            std::vector<std::tuple<unsigned int, unsigned int>> result{};
            for (auto it = pimpl->internalPotentialO2Registry.begin();
                 it != pimpl->internalPotentialO2Registry.end(); ++it) {
                result.push_back(std::make_tuple(it->first.t1, it->first.t2));
            }
            return result;
        }

        const boost::uuids::uuid &
        KernelContext::registerOrder1Potential(potentials::PotentialOrder1 const *const potential,
                                               const std::string &type) {
            auto typeId = pimpl->typeMapping[type];
            if (pimpl->internalPotentialO1Registry.find(typeId) == pimpl->internalPotentialO1Registry.end()) {
                pimpl->internalPotentialO1Registry.insert(
                        std::make_pair(typeId, std::vector<std::unique_ptr<potentials::PotentialOrder1>>()));
            }
            auto ptr = potential->replicate();
            ptr->configureForType(typeId);
            pimpl->internalPotentialO1Registry[typeId].push_back(std::unique_ptr<potentials::PotentialOrder1>(ptr));
            return pimpl->internalPotentialO1Registry[typeId].back()->getId();
        }

        std::vector<potentials::PotentialOrder1 *> KernelContext::getOrder1Potentials(const std::string &type) const {
            return getOrder1Potentials(pimpl->typeMapping[type]);
        }

        std::vector<potentials::PotentialOrder1 *> KernelContext::getOrder1Potentials(const unsigned int type) const {
            return readdy::utils::collections::getOrDefault(pimpl->potentialO1Registry, type, pimpl->defaultPotentialsO1);
        }

        std::unordered_set<unsigned int> KernelContext::getAllOrder1RegisteredPotentialTypes() const {
            std::unordered_set<unsigned int> result{};
            for (auto it = pimpl->internalPotentialO1Registry.begin();
                 it != pimpl->internalPotentialO1Registry.end(); ++it) {
                result.insert(it->first);
            }
            return result;
        }

        void KernelContext::deregisterPotential(const boost::uuids::uuid &potential) {
            const auto deleterO1 = [&potential](std::unique_ptr<potentials::PotentialOrder1> const &p) -> bool {
                return potential == p->getId();
            };
            const auto deleterO2 = [&potential](std::unique_ptr<potentials::PotentialOrder2> const &p) -> bool {
                return potential == p->getId();
            };
            for (auto it = pimpl->internalPotentialO1Registry.begin();
                 it != pimpl->internalPotentialO1Registry.end(); ++it) {
                it->second.erase(std::remove_if(it->second.begin(), it->second.end(), deleterO1), it->second.end());
            }
            for (auto it = pimpl->internalPotentialO2Registry.begin();
                 it != pimpl->internalPotentialO2Registry.end(); ++it) {
                it->second.erase(std::remove_if(it->second.begin(), it->second.end(), deleterO2), it->second.end());
            }
        }

        const std::vector<reactions::Reaction<1> *> &KernelContext::getOrder1Reactions(const std::string &type) const {
            return getOrder1Reactions(pimpl->typeMapping[type]);
        }

        const std::vector<reactions::Reaction<1> *> &KernelContext::getOrder1Reactions(const unsigned int type) const {
            return readdy::utils::collections::getOrDefault(*reactionOneEductRegistry, type, pimpl->defaultReactionsO1);
        }

        const std::vector<reactions::Reaction<2> *> &
        KernelContext::getOrder2Reactions(const std::string &type1, const std::string &type2) const {
            return getOrder2Reactions(pimpl->typeMapping[type1], pimpl->typeMapping[type2]);
        }

        const std::vector<reactions::Reaction<2> *> &
        KernelContext::getOrder2Reactions(const unsigned int type1, const unsigned int type2) const {
            return readdy::utils::collections::getOrDefault(*reactionTwoEductsRegistry, {type1, type2}, pimpl->defaultReactionsO2);
        }

        const std::vector<const reactions::Reaction<1> *> KernelContext::getAllOrder1Reactions() const {
            auto result = std::vector<const reactions::Reaction<1> *>();
            for (const auto &mapEntry : *internalReactionOneEductRegistry) {
                for (const auto &reaction : mapEntry.second) {
                    result.push_back(reaction.get());
                }
            }
            return result;
        }

        const reactions::Reaction<1> *const KernelContext::getReactionOrder1WithName(const std::string &name) const {
            for (const auto &mapEntry : *internalReactionOneEductRegistry) {
                for (const auto &reaction : mapEntry.second) {
                    if (reaction->getName() == name) return reaction.get();
                }
            }

            return nullptr;
        }

        const std::vector<const reactions::Reaction<2> *> KernelContext::getAllOrder2Reactions() const {
            auto result = std::vector<const reactions::Reaction<2> *>();
            for (const auto &mapEntry : *internalReactionTwoEductsRegistry) {
                for (const auto &reaction : mapEntry.second) {
                    result.push_back(reaction.get());
                }
            }
            return result;
        }

        const reactions::Reaction<2> *const KernelContext::getReactionOrder2WithName(const std::string &name) const {
            for (const auto &mapEntry : *internalReactionTwoEductsRegistry) {
                for (const auto &reaction : mapEntry.second) {
                    if (reaction->getName() == name) return reaction.get();
                }
            }
            return nullptr;
        }

        const std::function<void(Vec3 &)> &KernelContext::getFixPositionFun() const {
            return pimpl->fixPositionFun;
        }

        const std::function<double(const Vec3 &, const Vec3 &)> &KernelContext::getDistSquaredFun() const {
            return pimpl->distFun;
        }

        const std::function<Vec3(const Vec3 &, const Vec3 &)> &KernelContext::getShortestDifferenceFun() const {
            return pimpl->diffFun;
        }

        void KernelContext::configure() {
            pimpl->potentialO1Registry.clear();
            pimpl->potentialO2Registry.clear();
            reactionOneEductRegistry->clear();
            reactionTwoEductsRegistry->clear();
            for(auto typeIt : pimpl->typeMapping) {

            }
            for (auto &&e : pimpl->internalPotentialO1Registry) {
                pimpl->potentialO1Registry[e.first] = std::vector<readdy::model::potentials::PotentialOrder1 *>();
                pimpl->potentialO1Registry[e.first].reserve(e.second.size());
                for (auto &&pot : e.second) {
                    pot->configureForType(e.first);
                    pimpl->potentialO1Registry[e.first].push_back(pot.get());
                }
            }
            for (auto &&e : pimpl->internalPotentialO2Registry) {
                pimpl->potentialO2Registry[e.first] = std::vector<readdy::model::potentials::PotentialOrder2 *>();
                pimpl->potentialO2Registry[e.first].reserve(e.second.size());
                for (auto &&pot : e.second) {
                    pot->configureForTypes(e.first.t1, e.first.t2);
                    pimpl->potentialO2Registry[e.first].push_back(pot.get());
                }
            }
            for (auto &&e : *internalReactionOneEductRegistry) {
                (*reactionOneEductRegistry)[e.first] = std::vector<reactions::Reaction<1> *>();
                (*reactionOneEductRegistry)[e.first].reserve(e.second.size());
                std::for_each(e.second.begin(), e.second.end(),
                              [this, &e](const std::unique_ptr<reactions::Reaction<1>> &ptr) {
                                  (*reactionOneEductRegistry)[e.first].push_back(ptr.get());
                              }
                );
            }
            for (auto &&e : *internalReactionTwoEductsRegistry) {
                (*reactionTwoEductsRegistry)[e.first] = std::vector<reactions::Reaction<2> *>();
                (*reactionTwoEductsRegistry)[e.first].reserve(e.second.size());
                std::for_each(e.second.begin(), e.second.end(),
                              [this, &e](const std::unique_ptr<reactions::Reaction<2>> &ptr) {
                                  (*reactionTwoEductsRegistry)[e.first].push_back(ptr.get());
                              }
                );
            }
        }

        std::vector<unsigned int> KernelContext::getAllRegisteredParticleTypes() const {
            std::vector<unsigned int> v;
            for (auto &&entry : pimpl->typeMapping) {
                v.push_back(entry.second);
            }
            return v;
        }

        std::string KernelContext::getParticleName(unsigned int id) const {
            for (auto &&e : pimpl->typeMapping) {
                if (e.second == id) return e.first;
            }
            return "";
        }

        const std::unordered_map<unsigned int, std::vector<potentials::PotentialOrder1 *>>
        KernelContext::getAllOrder1Potentials() const {
            return pimpl->potentialO1Registry;
        }

        const std::unordered_map<readdy::util::ParticleTypePair, std::vector<potentials::PotentialOrder2 *>, readdy::util::ParticleTypePairHasher>
        KernelContext::getAllOrder2Potentials() const {
            return pimpl->potentialO2Registry;
        }

        const KernelContext::rdy_type_mapping &KernelContext::getTypeMapping() const {
            return pimpl->typeMapping;
        }

        const std::vector<potentials::PotentialOrder2 *> KernelContext::getVectorAllOrder2Potentials() const {
            std::vector<potentials::PotentialOrder2*> result;
            for (auto &&e : getAllOrder2RegisteredPotentialTypes()) {
                for (auto &&p : getOrder2Potentials(std::get<0>(e), std::get<1>(e))) {
                    result.push_back(p);
                }
            }
            return result;

        KernelContext &KernelContext::operator=(KernelContext &&rhs) = default;

        KernelContext::KernelContext(KernelContext &&rhs) = default;

        KernelContext::~KernelContext() = default;
    }
}






