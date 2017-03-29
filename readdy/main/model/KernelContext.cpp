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

KernelContext::KernelContext()
        : pimpl(std::make_unique<KernelContext::Impl>()),
          potentialRegistry_(std::cref(particleTypeRegistry_)),
          reactionRegistry_(std::cref(particleTypeRegistry_)) {}

std::array<double, 3> &KernelContext::getBoxSize() const {
    return pimpl->box_size;
}

const std::array<bool, 3> &KernelContext::getPeriodicBoundary() const {
    return pimpl->periodic_boundary;
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
    potentialRegistry_.configure();
    reactionRegistry_.configure();

    /**
     * Info output
     */
    if (debugOutput) {

        log::debug("Configured kernel context with: ");
        log::debug("--------------------------------");
        log::debug(" - kBT = {}", getKBT());
        log::debug(" - periodic b.c. = ({}, {}, {})", getPeriodicBoundary()[0], getPeriodicBoundary()[1],
                   getPeriodicBoundary()[2]);
        log::debug(" - box size = ({}, {}, {})", getBoxSize()[0], getBoxSize()[1], getBoxSize()[2]);

        potentialRegistry_.debug_output();
        reactionRegistry_.debug_output();
    }

}

std::tuple<readdy::model::Vec3, readdy::model::Vec3> KernelContext::getBoxBoundingVertices() const {
    const auto &boxSize = getBoxSize();
    readdy::model::Vec3 lowerLeft{-0.5 * boxSize[0], -0.5 * boxSize[1], -0.5 * boxSize[2]};
    readdy::model::Vec3 upperRight = lowerLeft + readdy::model::Vec3(boxSize);
    return std::make_tuple(std::move(lowerLeft), std::move(upperRight));
}

const KernelContext::compartment_registry &KernelContext::getCompartments() const {
    return *compartmentRegistry;
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

api::PotentialConfiguration &KernelContext::topologyPotentialConfiguration() {
    return potentialConfiguration_;
}

const api::PotentialConfiguration &KernelContext::topologyPotentialConfiguration() const {
    return potentialConfiguration_;
}

void KernelContext::configureTopologyBondPotential(const std::string &type1, const std::string &type2,
                                                   const api::Bond &bond) {
    potentialConfiguration_.pairPotentials[std::make_tuple(particleTypeRegistry_.getParticleTypeID(type1),
                                                           particleTypeRegistry_.getParticleTypeID(type2))].push_back(
            bond);
}

void KernelContext::configureTopologyAnglePotential(const std::string &type1, const std::string &type2,
                                                    const std::string &type3, const api::Angle &angle) {
    potentialConfiguration_.anglePotentials[std::make_tuple(particleTypeRegistry_.getParticleTypeID(type1),
                                                            particleTypeRegistry_.getParticleTypeID(type2),
                                                            particleTypeRegistry_.getParticleTypeID(type3))].push_back(
            angle);
}

void KernelContext::configureTopologyTorsionPotential(const std::string &type1, const std::string &type2,
                                                      const std::string &type3, const std::string &type4,
                                                      const api::TorsionAngle &torsionAngle) {
    potentialConfiguration_.torsionPotentials[std::make_tuple(particleTypeRegistry_.getParticleTypeID(type1),
                                                              particleTypeRegistry_.getParticleTypeID(type2),
                                                              particleTypeRegistry_.getParticleTypeID(type3),
                                                              particleTypeRegistry_.getParticleTypeID(
                                                                      type4))].push_back(torsionAngle);
}

reactions::ReactionRegistry &KernelContext::reactionRegistry() {
    return reactionRegistry_;
}

const reactions::ReactionRegistry &KernelContext::reactionRegistry() const {
    return reactionRegistry_;
}

ParticleTypeRegistry &KernelContext::particleTypeRegistry() {
    return particleTypeRegistry_;
}

const ParticleTypeRegistry &KernelContext::particleTypeRegistry() const {
    return particleTypeRegistry_;
}

const potentials::PotentialRegistry &KernelContext::potentialRegistry() const {
    return potentialRegistry_;
}

potentials::PotentialRegistry &KernelContext::potentialRegistry() {
    return potentialRegistry_;
}

KernelContext::~KernelContext() = default;


}
}






