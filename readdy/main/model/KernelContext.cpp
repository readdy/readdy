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
    unsigned int typeCounter;
    std::unordered_map<std::string, unsigned int> typeMapping;
    double kBT = 1;
    std::array<double, 3> box_size{{1, 1, 1}};
    std::array<bool, 3> periodic_boundary{{true, true, true}};
    std::unordered_map<unsigned int, double> diffusionConstants{};
    std::unordered_map<unsigned int, double> particleRadii{};

    double timeStep;

    std::function<void(Vec3 &)> fixPositionFun = [](
            Vec3 &vec) -> void { readdy::model::fixPosition<true, true, true>(vec, 1., 1., 1.); };
    std::function<Vec3(const Vec3 &, const Vec3 &)> diffFun = [](const Vec3 &lhs, const Vec3 &rhs) -> Vec3 {
        return readdy::model::shortestDifference<true, true, true>(lhs, rhs, 1., 1., 1.);
    };
    std::function<double(const Vec3 &, const Vec3 &)> distFun = [&](const Vec3 &lhs,
                                                                    const Vec3 &rhs) -> double {
        const auto dv = diffFun(lhs, rhs);
        return dv * dv;
    };

    std::vector<reactions::Reaction<1> *> defaultReactionsO1{};
    std::vector<reactions::Reaction<2> *> defaultReactionsO2{};
    std::vector<potentials::PotentialOrder1 *> defaultPotentialsO1{};
    std::vector<potentials::PotentialOrder2 *> defaultPotentialsO2{};

    ~Impl() = default;

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

// todo respect const correctness, dont create new entries, thx
unsigned int KernelContext::getParticleTypeID(const std::string &name) const {
    return pimpl->typeMapping[name];
}

double KernelContext::getDiffusionConstant(unsigned int particleType) const {
    return pimpl->diffusionConstants[particleType];
}

double KernelContext::getParticleRadius(const std::string &type) const {
    return getParticleRadius(pimpl->typeMapping[type]);
}

double KernelContext::getParticleRadius(const unsigned int type) const {
    if (pimpl->particleRadii.find(type) == pimpl->particleRadii.end()) {
        log::console()->warn("No particle radius was set for the particle type id {}, setting r=1", type);
        pimpl->particleRadii[type] = 1;
    }
    return pimpl->particleRadii[type];
}

void KernelContext::setParticleRadius(const std::string &particleType, const double r) {
    pimpl->particleRadii[getOrCreateTypeId(particleType)] = r;
}

unsigned int KernelContext::getOrCreateTypeId(const std::string &particleType) {
    unsigned int t_id;
    if (pimpl->typeMapping.find(particleType) != pimpl->typeMapping.end()) {
        t_id = pimpl->typeMapping[particleType];
    } else {
        t_id = ++(pimpl->typeCounter);
        pimpl->typeMapping.emplace(particleType, t_id);
    }
    return t_id;
}

const std::vector<potentials::PotentialOrder2 *> &
KernelContext::getOrder2Potentials(const std::string &type1, const std::string &type2) const {
    return getOrder2Potentials(pimpl->typeMapping[type1], pimpl->typeMapping[type2]);
}

const std::vector<potentials::PotentialOrder2 *> &
KernelContext::getOrder2Potentials(const unsigned int type1, const unsigned int type2) const {
    return readdy::util::collections::getOrDefault(*potentialO2Registry, {type1, type2},
                                                   pimpl->defaultPotentialsO2);
}

std::vector<std::tuple<unsigned int, unsigned int>>
KernelContext::getAllOrder2RegisteredPotentialTypes() const {
    std::vector<std::tuple<unsigned int, unsigned int>> result{};
    for (auto it = potentialO2Registry->begin();
         it != potentialO2Registry->end(); ++it) {
        result.push_back(std::make_tuple(it->first.t1, it->first.t2));
    }
    return result;
}

std::vector<potentials::PotentialOrder1 *> KernelContext::getOrder1Potentials(const std::string &type) const {
    return getOrder1Potentials(pimpl->typeMapping[type]);
}

std::vector<potentials::PotentialOrder1 *> KernelContext::getOrder1Potentials(const unsigned int type) const {
    return readdy::util::collections::getOrDefault(*potentialO1Registry, type, pimpl->defaultPotentialsO1);
}

std::unordered_set<unsigned int> KernelContext::getAllOrder1RegisteredPotentialTypes() const {
    std::unordered_set<unsigned int> result{};
    for (auto it = potentialO1RegistryInternal->begin();
         it != potentialO1RegistryInternal->end(); ++it) {
        result.insert(it->first);
    }
    return result;
}

void KernelContext::deregisterPotential(const short potential) {
    for (auto it = potentialO1RegistryInternal->begin(); it != potentialO1RegistryInternal->end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [&potential](const std::unique_ptr<potentials::PotentialOrder1> &p) -> bool {
                                            return potential == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO2RegistryInternal->begin(); it != potentialO2RegistryInternal->end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [&potential](const std::unique_ptr<potentials::PotentialOrder2> &p) -> bool {
                                            return potential == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO1RegistryExternal->begin(); it != potentialO1RegistryExternal->end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [&potential](potentials::PotentialOrder1* p) -> bool {
                                            return potential == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO2RegistryExternal->begin(); it != potentialO2RegistryExternal->end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [&potential](potentials::PotentialOrder2* p) -> bool {
                                            return potential == p->getId();
                                        }
        ), it->second.end());
    }
}

const std::vector<reactions::Reaction<1> *> &KernelContext::getOrder1Reactions(const std::string &type) const {
    return getOrder1Reactions(pimpl->typeMapping[type]);
}

const std::vector<reactions::Reaction<1> *> &KernelContext::getOrder1Reactions(const unsigned int type) const {
    return readdy::util::collections::getOrDefault(*reactionOneEductRegistry, type, pimpl->defaultReactionsO1);
}

const std::vector<reactions::Reaction<2> *> &
KernelContext::getOrder2Reactions(const std::string &type1, const std::string &type2) const {
    return getOrder2Reactions(pimpl->typeMapping[type1], pimpl->typeMapping[type2]);
}

const std::vector<reactions::Reaction<2> *> &
KernelContext::getOrder2Reactions(const unsigned int type1, const unsigned int type2) const {
    return readdy::util::collections::getOrDefault(*reactionTwoEductsRegistry, {type1, type2},
                                                   pimpl->defaultReactionsO2);
}

const std::vector<const reactions::Reaction<1> *> KernelContext::getAllOrder1Reactions() const {
    auto result = std::vector<const reactions::Reaction<1> *>();
    for (const auto &mapEntry : *reactionOneEductRegistry) {
        for (const auto &reaction : mapEntry.second) {
            result.push_back(reaction);
        }
    }
    return result;
}

const reactions::Reaction<1> *const KernelContext::getReactionOrder1WithName(const std::string &name) const {
    for (const auto &mapEntry : *reactionOneEductRegistry) {
        for (const auto &reaction : mapEntry.second) {
            if (reaction->getName() == name) return reaction;
        }
    }

    return nullptr;
}

const std::vector<const reactions::Reaction<2> *> KernelContext::getAllOrder2Reactions() const {
    auto result = std::vector<const reactions::Reaction<2> *>();
    for (const auto &mapEntry : *reactionTwoEductsRegistry) {
        for (const auto &reaction : mapEntry.second) {
            result.push_back(reaction);
        }
    }
    return result;
}

const reactions::Reaction<2> *const KernelContext::getReactionOrder2WithName(const std::string &name) const {
    for (const auto &mapEntry : *reactionTwoEductsRegistry) {
        for (const auto &reaction : mapEntry.second) {
            if (reaction->getName() == name) return reaction;
        }
    }
    return nullptr;
}

const KernelContext::fix_pos_fun &KernelContext::getFixPositionFun() const {
    return pimpl->fixPositionFun;
}

const KernelContext::dist_squared_fun &KernelContext::getDistSquaredFun() const {
    return pimpl->distFun;
}

const KernelContext::shortest_dist_fun &KernelContext::getShortestDifferenceFun() const {
    return pimpl->diffFun;
}

void KernelContext::configure(bool debugOutput) {
    namespace coll = readdy::util::collections;
    using pair = util::ParticleTypePair;
    using pot1 = potentials::PotentialOrder1;
    using pot1_ptr = std::unique_ptr<potentials::PotentialOrder1>;
    using pot2_ptr = std::unique_ptr<potentials::PotentialOrder2>;
    using pot2 = potentials::PotentialOrder2;
    using reaction1ptr = std::unique_ptr<reactions::Reaction<1>>;
    using reaction2ptr = std::unique_ptr<reactions::Reaction<2>>;

    potentialO1Registry->clear();
    potentialO2Registry->clear();
    reactionOneEductRegistry->clear();
    reactionTwoEductsRegistry->clear();

    coll::for_each_value(*potentialO1RegistryInternal, [&](const unsigned int type, const pot1_ptr& ptr) {
        ptr->configureForType(type); (*potentialO1Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(*potentialO2RegistryInternal, [&](const pair& type, const pot2_ptr& ptr) {
        ptr->configureForTypes(type.t1, type.t2); (*potentialO2Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(*potentialO1RegistryExternal, [&](const unsigned int type, pot1* ptr) {
        ptr->configureForType(type); (*potentialO1Registry)[type].push_back(ptr);
    });
    coll::for_each_value(*potentialO2RegistryExternal, [&](const pair& type, pot2* ptr) {
        ptr->configureForTypes(type.t1, type.t2); (*potentialO2Registry)[type].push_back(ptr);
    });
    coll::for_each_value(*reactionOneEductRegistryInternal, [&](const unsigned int type, const reaction1ptr& ptr) {
        (*reactionOneEductRegistry)[type].push_back(ptr.get());
    });
    coll::for_each_value(*reactionTwoEductsRegistryInternal, [&](const pair& type, const reaction2ptr& r) {
        (*reactionTwoEductsRegistry)[type].push_back(r.get());
    });
    coll::for_each_value(*reactionOneEductRegistryExternal, [&](const unsigned int type, reactions::Reaction<1>* ptr) {
        (*reactionOneEductRegistry)[type].push_back(ptr);
    });
    coll::for_each_value(*reactionTwoEductsRegistryExternal, [&](const pair& type, reactions::Reaction<2>* r) {
        (*reactionTwoEductsRegistry)[type].push_back(r);
    });

    /**
     * Info output
     */
    if(debugOutput) {
        auto find_pot_name = [this](unsigned int type) -> const std::string {
            for (auto &&t : pimpl->typeMapping) {
                if (t.second == type) return t.first;
            }
            return "";
        };
        log::console()->debug("Configured kernel context with: ");
        log::console()->debug("--------------------------------");
        log::console()->debug(" - kBT = {}", getKBT());
        log::console()->debug(" - tau = {}", getTimeStep());
        log::console()->debug(" - periodic b.c. = ({}, {}, {})", getPeriodicBoundary()[0], getPeriodicBoundary()[1],
                              getPeriodicBoundary()[2]);
        log::console()->debug(" - box size = ({}, {}, {})", getBoxSize()[0], getBoxSize()[1], getBoxSize()[2]);

        if (!getAllOrder1Potentials().empty()) {
            log::console()->debug(" - potentials of order 1:");
            for (auto types : getAllOrder1Potentials()) {
                log::console()->debug("     * for type {}", find_pot_name(types.first));
                for (auto pot : types.second) {
                    log::console()->debug("         * {}", pot->getName());
                }
            }
        }
        if (!getAllOrder2Potentials().empty()) {
            log::console()->debug(" - potentials of order 2:");
            for (auto types : getAllOrder2Potentials()) {
                log::console()->debug("     * for types {} and {}", find_pot_name(types.first.t1),
                                      find_pot_name(types.first.t2));
                for (auto pot : types.second) {
                    log::console()->debug("         * {}", pot->getName());
                }
            }
        }
        if (!getAllOrder1Reactions().empty()) {
            log::console()->debug(" - reactions of order 1:");
            for (auto reaction : getAllOrder1Reactions()) {
                log::console()->debug("     * reaction {}", *reaction);
            }
        }
        if (!getAllOrder2Reactions().empty()) {
            log::console()->debug(" - reactions of order 2:");
            for (auto reaction : getAllOrder2Reactions()) {
                log::console()->debug("     * reaction {}", *reaction);
            }
        }
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
    return *potentialO1Registry;
}

const std::unordered_map<readdy::util::ParticleTypePair, std::vector<potentials::PotentialOrder2 *>, readdy::util::ParticleTypePairHasher>
KernelContext::getAllOrder2Potentials() const {
    return *potentialO2Registry;
}

const KernelContext::rdy_type_mapping &KernelContext::getTypeMapping() const {
    return pimpl->typeMapping;
}

const std::vector<potentials::PotentialOrder2 *> KernelContext::getVectorAllOrder2Potentials() const {
    std::vector<potentials::PotentialOrder2 *> result;
    for (auto &&e : getAllOrder2RegisteredPotentialTypes()) {
        for (auto &&p : getOrder2Potentials(std::get<0>(e), std::get<1>(e))) {
            result.push_back(p);
        }
    }
    return result;
}

std::tuple<readdy::model::Vec3, readdy::model::Vec3> KernelContext::getBoxBoundingVertices() const {
    const auto &boxSize = getBoxSize();
    readdy::model::Vec3 lowerLeft{-0.5 * boxSize[0], -0.5 * boxSize[1], -0.5 * boxSize[2]};
    readdy::model::Vec3 upperRight = lowerLeft + readdy::model::Vec3(boxSize);
    return std::make_tuple(std::move(lowerLeft), std::move(upperRight));
}

KernelContext &KernelContext::operator=(KernelContext &&rhs) = default;

KernelContext::KernelContext(KernelContext &&rhs) = default;

KernelContext::~KernelContext() = default;
}
}






