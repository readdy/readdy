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
 * @file KernelContext.cpp
 * @brief Implementation file of the KernelContext.
 * @author clonker
 * @author chrisfroe
 * @date 18.04.16
 * @todo make proper reference to KernelContext.h, is kBT really indepdendent of t?
 */

#include <readdy/model/KernelContext.h>
#include <readdy/common/Utils.h>
#include <readdy/model/_internal/Util.h>

namespace readdy {
namespace model {

using particle_t = readdy::model::Particle;

struct KernelContext::Impl {

    particle_t::type_type typeCounter = 0;
    std::unordered_map<std::string, particle_t::type_type> typeMapping;
    std::unordered_map<particle_t::type_type, ParticleTypeInfo> particleInfo;

    double kBT = 1;
    std::array<double, 3> box_size{{1, 1, 1}};
    std::array<bool, 3> periodic_boundary{{true, true, true}};

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

void KernelContext::registerParticleType(const std::string &name, const double diffusionConst, const double radius,
                                         const readdy::model::Particle::flavor_t flavor) {
    particle_t::type_type t_id = (pimpl->typeCounter)++;
    pimpl->typeMapping.emplace(name, t_id);
    pimpl->particleInfo.emplace(std::make_pair(t_id, ParticleTypeInfo{name, diffusionConst, radius, flavor, t_id}));
}

double KernelContext::getDiffusionConstant(const std::string &particleType) const {
    return getDiffusionConstant(getParticleTypeID(particleType));
}

double KernelContext::getDiffusionConstant(particle_t::type_type particleType) const {
    return pimpl->particleInfo.at(particleType).diffusionConstant;
}

particle_t::type_type KernelContext::getParticleTypeID(const std::string &name) const {
    return pimpl->typeMapping.at(name);
}

double KernelContext::getParticleRadius(const std::string &particleType) const {
    return getParticleRadius(getParticleTypeID(particleType));
}

double KernelContext::getParticleRadius(const particle_t::type_type type) const {
    return pimpl->particleInfo.at(type).radius;
}

const std::vector<potentials::PotentialOrder2 *> &
KernelContext::getOrder2Potentials(const std::string &type1, const std::string &type2) const {
    return getOrder2Potentials(pimpl->typeMapping[type1], pimpl->typeMapping[type2]);
}

const std::vector<potentials::PotentialOrder2 *> &
KernelContext::getOrder2Potentials(const particle_t::type_type type1, const particle_t::type_type type2) const {
    return readdy::util::collections::getOrDefault(potentialO2Registry, std::make_tuple(type1, type2),
                                                   pimpl->defaultPotentialsO2);
}

std::vector<std::tuple<particle_t::type_type, particle_t::type_type>>
KernelContext::getAllOrder2RegisteredPotentialTypes() const {
    std::vector<std::tuple<particle_t::type_type, particle_t::type_type>> result{};
    for (auto it = potentialO2Registry.begin();
         it != potentialO2Registry.end(); ++it) {
        result.push_back(it->first);
    }
    return result;
}

std::vector<potentials::PotentialOrder1 *> KernelContext::getOrder1Potentials(const std::string &type) const {
    return getOrder1Potentials(pimpl->typeMapping[type]);
}

std::vector<potentials::PotentialOrder1 *> KernelContext::getOrder1Potentials(const particle_t::type_type type) const {
    return readdy::util::collections::getOrDefault(potentialO1Registry, type, pimpl->defaultPotentialsO1);
}

std::unordered_set<particle_t::type_type> KernelContext::getAllOrder1RegisteredPotentialTypes() const {
    std::unordered_set<particle_t::type_type> result{};
    for (auto it = potentialO1RegistryInternal.begin();
         it != potentialO1RegistryInternal.end(); ++it) {
        result.insert(it->first);
    }
    return result;
}

void KernelContext::deregisterPotential(const short potential) {
    for (auto it = potentialO1RegistryInternal.begin(); it != potentialO1RegistryInternal.end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [&potential](const std::unique_ptr<potentials::PotentialOrder1> &p) -> bool {
                                            return potential == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO2RegistryInternal.begin(); it != potentialO2RegistryInternal.end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [&potential](const std::unique_ptr<potentials::PotentialOrder2> &p) -> bool {
                                            return potential == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO1RegistryExternal.begin(); it != potentialO1RegistryExternal.end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [&potential](potentials::PotentialOrder1 *p) -> bool {
                                            return potential == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO2RegistryExternal.begin(); it != potentialO2RegistryExternal.end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [&potential](potentials::PotentialOrder2 *p) -> bool {
                                            return potential == p->getId();
                                        }
        ), it->second.end());
    }
}

const std::vector<reactions::Reaction<1> *> &KernelContext::getOrder1Reactions(const std::string &type) const {
    return getOrder1Reactions(pimpl->typeMapping[type]);
}

const std::vector<reactions::Reaction<1> *> &KernelContext::getOrder1Reactions(const particle_t::type_type type) const {
    return readdy::util::collections::getOrDefault(reactionOneEductRegistry, type, pimpl->defaultReactionsO1);
}

const std::vector<reactions::Reaction<2> *> &
KernelContext::getOrder2Reactions(const std::string &type1, const std::string &type2) const {
    return getOrder2Reactions(pimpl->typeMapping[type1], pimpl->typeMapping[type2]);
}

const std::vector<reactions::Reaction<2> *> &
KernelContext::getOrder2Reactions(const particle_t::type_type type1, const particle_t::type_type type2) const {
    return readdy::util::collections::getOrDefault(reactionTwoEductsRegistry, std::make_tuple(type1, type2),
                                                   pimpl->defaultReactionsO2);
}

const std::vector<const reactions::Reaction<1> *> KernelContext::getAllOrder1Reactions() const {
    auto result = std::vector<const reactions::Reaction<1> *>();
    for (const auto &mapEntry : reactionOneEductRegistry) {
        for (const auto reaction : mapEntry.second) {
            result.push_back(reaction);
        }
    }
    return result;
}

const reactions::Reaction<1> *const KernelContext::getReactionOrder1WithName(const std::string &name) const {
    for (const auto &mapEntry : reactionOneEductRegistry) {
        for (const auto &reaction : mapEntry.second) {
            if (reaction->getName() == name) return reaction;
        }
    }

    return nullptr;
}

const std::vector<const reactions::Reaction<2> *> KernelContext::getAllOrder2Reactions() const {
    auto result = std::vector<const reactions::Reaction<2> *>();
    for (const auto &mapEntry : reactionTwoEductsRegistry) {
        for (const auto reaction : mapEntry.second) {
            result.push_back(reaction);
        }
    }
    return result;
}

const reactions::Reaction<2> *const KernelContext::getReactionOrder2WithName(const std::string &name) const {
    for (const auto &mapEntry : reactionTwoEductsRegistry) {
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
    using pair = util::particle_type_pair;
    using pot1 = potentials::PotentialOrder1;
    using pot1_ptr = std::unique_ptr<potentials::PotentialOrder1>;
    using pot2_ptr = std::unique_ptr<potentials::PotentialOrder2>;
    using pot2 = potentials::PotentialOrder2;
    using reaction1ptr = std::unique_ptr<reactions::Reaction<1>>;
    using reaction2ptr = std::unique_ptr<reactions::Reaction<2>>;

    potentialO1Registry.clear();
    potentialO2Registry.clear();
    reactionOneEductRegistry.clear();
    reactionTwoEductsRegistry.clear();

    coll::for_each_value(potentialO1RegistryInternal, [&](const particle_t::type_type type, const pot1_ptr &ptr) {
        ptr->configureForType(this, type);
        (potentialO1Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(potentialO2RegistryInternal, [&](const pair &type, const pot2_ptr &ptr) {
        ptr->configureForTypes(this, std::get<0>(type), std::get<1>(type));
        (potentialO2Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(potentialO1RegistryExternal, [&](const particle_t::type_type type, pot1 *ptr) {
        ptr->configureForType(this, type);
        (potentialO1Registry)[type].push_back(ptr);
    });
    coll::for_each_value(potentialO2RegistryExternal, [&](const pair &type, pot2 *ptr) {
        ptr->configureForTypes(this, std::get<0>(type), std::get<1>(type));
        (potentialO2Registry)[type].push_back(ptr);
    });
    coll::for_each_value(reactionOneEductRegistryInternal,
                         [&](const particle_t::type_type type, const reaction1ptr &ptr) {
                             (reactionOneEductRegistry)[type].push_back(ptr.get());
                         });
    coll::for_each_value(reactionTwoEductsRegistryInternal, [&](const pair &type, const reaction2ptr &r) {
        (reactionTwoEductsRegistry)[type].push_back(r.get());
    });
    coll::for_each_value(reactionOneEductRegistryExternal,
                         [&](const particle_t::type_type type, reactions::Reaction<1> *ptr) {
                             (reactionOneEductRegistry)[type].push_back(ptr);
                         });
    coll::for_each_value(reactionTwoEductsRegistryExternal, [&](const pair &type, reactions::Reaction<2> *r) {
        (reactionTwoEductsRegistry)[type].push_back(r);
    });

    /**
     * Info output
     */
    if (debugOutput) {
        auto find_pot_name = [this](particle_t::type_type type) -> const std::string {
            for (auto &&t : pimpl->typeMapping) {
                if (t.second == type) return t.first;
            }
            return "";
        };
        log::debug("Configured kernel context with: ");
        log::debug("--------------------------------");
        log::debug(" - kBT = {}", getKBT());
        log::debug(" - periodic b.c. = ({}, {}, {})", getPeriodicBoundary()[0], getPeriodicBoundary()[1],
                   getPeriodicBoundary()[2]);
        log::debug(" - box size = ({}, {}, {})", getBoxSize()[0], getBoxSize()[1], getBoxSize()[2]);

        if (!getAllOrder1Potentials().empty()) {
            log::debug(" - potentials of order 1:");
            for (const auto &types : getAllOrder1Potentials()) {
                log::debug("     * for type {}", find_pot_name(types.first));
                for (auto pot : types.second) {
                    log::debug("         * {}", _internal::util::to_string(pot));
                }
            }
        }
        if (!getAllOrder2Potentials().empty()) {
            log::debug(" - potentials of order 2:");
            for (const auto &types : getAllOrder2Potentials()) {
                log::debug("     * for types {} and {}", find_pot_name(std::get<0>(types.first)),
                           find_pot_name(std::get<1>(types.first)));
                for (auto pot : types.second) {
                    log::debug("         * {}", _internal::util::to_string(pot));
                }
            }
        }
        if (!getAllOrder1Reactions().empty()) {
            log::debug(" - reactions of order 1:");
            for (const auto &reaction : getAllOrder1Reactions()) {
                log::debug("     * reaction {}", *reaction);
            }
        }
        if (!getAllOrder2Reactions().empty()) {
            log::debug(" - reactions of order 2:");
            for (const auto &reaction : getAllOrder2Reactions()) {
                log::debug("     * reaction {}", *reaction);
            }
        }
    }

}

std::vector<particle_t::type_type> KernelContext::getAllRegisteredParticleTypes() const {
    std::vector<particle_t::type_type> v;
    for (auto &&entry : pimpl->typeMapping) {
        v.push_back(entry.second);
    }
    return v;
}

std::string KernelContext::getParticleName(particle_t::type_type id) const {
    for (auto &&e : pimpl->typeMapping) {
        if (e.second == id) return e.first;
    }
    return "";
}

const std::unordered_map<particle_t::type_type, std::vector<potentials::PotentialOrder1 *>>
KernelContext::getAllOrder1Potentials() const {
    return potentialO1Registry;
}

const KernelContext::rdy_pot_2_registry KernelContext::getAllOrder2Potentials() const {
    return potentialO2Registry;
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

const KernelContext::rdy_reverse_type_mapping KernelContext::generateReverseTypeMapping() const {
    const auto &typeMapping = getTypeMapping();
    rdy_reverse_type_mapping reverseTypeMapping;
    auto it = typeMapping.cbegin();
    while (it != typeMapping.cend()) {
        auto key = (*it).first;
        auto value = (*it).second;
        reverseTypeMapping.emplace(std::make_pair(value, key));
        ++it;
    }
    return reverseTypeMapping;
}

const KernelContext::compartment_registry &KernelContext::getCompartments() const {
    return *compartmentRegistry;
}

const ParticleTypeInfo &KernelContext::getParticleTypeInfo(const std::string &name) const {
    return getParticleTypeInfo(pimpl->typeMapping.at(name));
}

const ParticleTypeInfo &KernelContext::getParticleTypeInfo(const particle_t::type_type type) const {
    return pimpl->particleInfo.at(type);
}

const bool &KernelContext::recordReactionsWithPositions() const {
    return recordReactionsWithPositions_;
}

bool &KernelContext::recordReactionsWithPositions() {
    return recordReactionsWithPositions_;
}

const bool &KernelContext::recordReactionCounts() const {
    return recordReactionCounts_;
}

bool &KernelContext::recordReactionCounts() {
    return recordReactionCounts_;
}

top::graph::PotentialConfiguration &KernelContext::topologyPotentialConfiguration() {
    return potentialConfiguration_;
}

const top::graph::PotentialConfiguration &KernelContext::topologyPotentialConfiguration() const {
    return potentialConfiguration_;
}

void KernelContext::configureTopologyBondPotential(const std::string &type1, const std::string &type2,
                                                   const top::graph::Bond &bond) {
    potentialConfiguration_.pairPotentials[std::make_tuple(getParticleTypeID(type1),
                                                           getParticleTypeID(type2))].push_back(bond);
}

void KernelContext::configureTopologyAnglePotential(const std::string &type1, const std::string &type2,
                                                    const std::string &type3, const top::graph::Angle &angle) {
    potentialConfiguration_.anglePotentials[std::make_tuple(getParticleTypeID(type1), getParticleTypeID(type2),
                                                            getParticleTypeID(type3))].push_back(angle);
}

void KernelContext::configureTopologyTorsionPotential(const std::string &type1, const std::string &type2,
                                                      const std::string &type3, const std::string &type4,
                                                      const top::graph::TorsionAngle &torsionAngle) {
    potentialConfiguration_.torsionPotentials[std::make_tuple(getParticleTypeID(type1), getParticleTypeID(type2),
                                                              getParticleTypeID(type3),
                                                              getParticleTypeID(type4))].push_back(torsionAngle);
}

KernelContext &KernelContext::operator=(KernelContext &&rhs) = default;

KernelContext::KernelContext(KernelContext &&rhs) = default;

KernelContext::~KernelContext() = default;

ParticleTypeInfo::ParticleTypeInfo(const std::string &name, const double diffusionConstant, const double radius,
                                   const Particle::flavor_t flavor, const Particle::type_type typeId)
        : name(name), diffusionConstant(diffusionConstant), radius(radius), flavor(flavor), typeId(typeId) {}
}
}






