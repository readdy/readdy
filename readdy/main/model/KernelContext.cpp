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

const scalar &KernelContext::kBT() const {
    return _kBT;
}

scalar &KernelContext::kBT() {
    return _kBT;
}

KernelContext::KernelContext()
        : _potentialRegistry(_particleTypeRegistry), _reactionRegistry(_particleTypeRegistry),
          _topologyRegistry(_particleTypeRegistry) {
    using namespace std::placeholders;
    _pbc = std::bind(&readdy::model::applyPBC<false, false, false>, _1, c_::one, c_::one, c_::one);
    _fixPositionFun = std::bind(&readdy::model::fixPosition<false, false, false>, _1, c_::one, c_::one, c_::one);
    _diffFun = std::bind(&readdy::model::shortestDifference<false, false, false>, _1, _2, c_::one, c_::one, c_::one);
    _distFun = [this](const Vec3 &v1, const Vec3 &v2) {
        const auto dv = _diffFun(v1, v2);
        return dv * dv;
    };
}

const KernelContext::fix_pos_fun &KernelContext::fixPositionFun() const {
    return _fixPositionFun;
}

const KernelContext::dist_squared_fun &KernelContext::distSquaredFun() const {
    return _distFun;
}

const KernelContext::shortest_dist_fun &KernelContext::shortestDifferenceFun() const {
    return _diffFun;
}

void KernelContext::configure(bool debugOutput) {
    updateFunctions();

    _particleTypeRegistry.configure();
    _potentialRegistry.configure();
    _reactionRegistry.configure();
    _topologyRegistry.configure();

    /**
     * Info output
     */
    if (debugOutput) {

        log::debug("Configured kernel context with: ");
        log::debug("--------------------------------");
        log::debug(" - kBT = {}", kBT());
        log::debug(" - periodic b.c. = ({}, {}, {})", periodicBoundaryConditions()[0], periodicBoundaryConditions()[1],
                   periodicBoundaryConditions()[2]);
        log::debug(" - box size = ({}, {}, {})", boxSize()[0], boxSize()[1], boxSize()[2]);

        _particleTypeRegistry.debug_output();
        _potentialRegistry.debug_output();
        _reactionRegistry.debug_output();
        _topologyRegistry.debug_output();
    }

}

std::tuple<Vec3, Vec3> KernelContext::getBoxBoundingVertices() const {
    const auto &boxSize = _box_size;
    Vec3 lowerLeft{static_cast<scalar>(-0.5) * boxSize[0],
                   static_cast<scalar>(-0.5) * boxSize[1],
                   static_cast<scalar>(-0.5) * boxSize[2]};
    auto upperRight = lowerLeft + Vec3(boxSize);
    return std::make_tuple(lowerLeft, upperRight);
}

const KernelContext::CompartmentRegistry &KernelContext::compartments() const {
    return _compartmentRegistry;
}


const bool &KernelContext::recordReactionsWithPositions() const {
    return _recordReactionsWithPositions;
}

bool &KernelContext::recordReactionsWithPositions() {
    return _recordReactionsWithPositions;
}

const bool &KernelContext::recordReactionCounts() const {
    return _recordReactionCounts;
}

bool &KernelContext::recordReactionCounts() {
    return _recordReactionCounts;
}

reactions::ReactionRegistry &KernelContext::reactions() {
    return _reactionRegistry;
}

const reactions::ReactionRegistry &KernelContext::reactions() const {
    return _reactionRegistry;
}

ParticleTypeRegistry &KernelContext::particle_types() {
    return _particleTypeRegistry;
}

const ParticleTypeRegistry &KernelContext::particle_types() const {
    return _particleTypeRegistry;
}

const potentials::PotentialRegistry &KernelContext::potentials() const {
    return _potentialRegistry;
}

potentials::PotentialRegistry &KernelContext::potentials() {
    return _potentialRegistry;
}

const KernelContext::pbc_fun &KernelContext::applyPBCFun() const {
    return _pbc;
}

const scalar KernelContext::calculateMaxCutoff() const {
    scalar max_cutoff{0};
    for (const auto &entry : potentials().potentials_order2()) {
        for (const auto &potential : entry.second) {
            max_cutoff = std::max(max_cutoff, potential->getCutoffRadius());
        }
    }
    for (const auto &entry : reactions().order2()) {
        for (const auto &reaction : entry.second) {
            max_cutoff = std::max(max_cutoff, reaction->getEductDistance());
        }
    }
    for (const auto &entry : _topologyRegistry.spatial_reaction_registry()) {
        for (const auto &reaction : entry.second) {
            max_cutoff = std::max(max_cutoff, reaction.radius());
        }
    }
    return max_cutoff;
}

top::TopologyRegistry &KernelContext::topology_registry() {
    return _topologyRegistry;
}

const top::TopologyRegistry &KernelContext::topology_registry() const {
    return _topologyRegistry;
}

void KernelContext::updateFunctions() {
    using namespace std::placeholders;
    const auto &box = _box_size;
    if (_periodic_boundary[0]) {
        if (_periodic_boundary[1]) {
            if (_periodic_boundary[2]) {
                _pbc = std::bind(&applyPBC<true, true, true>, _1, box[0], box[1], box[2]);
                _fixPositionFun = std::bind(&fixPosition<true, true, true>, _1, box[0], box[1], box[2]);
                _diffFun = std::bind(&shortestDifference<true, true, true>, _1, _2, box[0], box[1], box[2]);
            } else {
                _pbc = std::bind(&applyPBC<true, true, false>, _1, box[0], box[1], box[2]);
                _fixPositionFun = std::bind(&fixPosition<true, true, false>, _1, box[0], box[1], box[2]);
                _diffFun = std::bind(&shortestDifference<true, true, false>, _1, _2, box[0], box[1], box[2]);
            }
        } else {
            if (_periodic_boundary[2]) {
                _pbc = std::bind(&applyPBC<true, false, true>, _1, box[0], box[1], box[2]);
                _fixPositionFun = std::bind(&fixPosition<true, false, true>, _1, box[0], box[1], box[2]);
                _diffFun = std::bind(&shortestDifference<true, false, true>, _1, _2, box[0], box[1], box[2]);
            } else {
                _pbc = std::bind(&applyPBC<true, false, false>, _1, box[0], box[1], box[2]);
                _fixPositionFun = std::bind(&fixPosition<true, false, false>, _1, box[0], box[1], box[2]);
                _diffFun = std::bind(&shortestDifference<true, false, false>, _1, _2, box[0], box[1], box[2]);
            }
        }
    } else {
        if (_periodic_boundary[1]) {
            if (_periodic_boundary[2]) {
                _pbc = std::bind(&applyPBC<false, true, true>, _1, box[0], box[1], box[2]);
                _fixPositionFun = std::bind(&fixPosition<false, true, true>, _1, box[0], box[1], box[2]);
                _diffFun = std::bind(&shortestDifference<false, true, true>, _1, _2, box[0], box[1], box[2]);
            } else {
                _pbc = std::bind(&applyPBC<false, true, false>, _1, box[0], box[1], box[2]);
                _fixPositionFun = std::bind(&fixPosition<false, true, false>, _1, box[0], box[1], box[2]);
                _diffFun = std::bind(&shortestDifference<false, true, false>, _1, _2, box[0], box[1], box[2]);
            }
        } else {
            if (_periodic_boundary[2]) {
                _pbc = std::bind(&applyPBC<false, false, true>, _1, box[0], box[1], box[2]);
                _fixPositionFun = std::bind(&fixPosition<false, false, true>, _1, box[0], box[1], box[2]);
                _diffFun = std::bind(&shortestDifference<false, false, true>, _1, _2, box[0], box[1], box[2]);
            } else {
                _pbc = std::bind(&applyPBC<false, false, false>, _1, box[0], box[1], box[2]);
                _fixPositionFun = std::bind(&fixPosition<false, false, false>, _1, box[0], box[1], box[2]);
                _diffFun = std::bind(&shortestDifference<false, false, false>, _1, _2, box[0], box[1], box[2]);
            }
        }
    }
}

const KernelContext::BoxSize &KernelContext::boxSize() const {
    return _box_size;
}

KernelContext::BoxSize &KernelContext::boxSize() {
    return _box_size;
}

const KernelContext::PeriodicBoundaryConditions &KernelContext::periodicBoundaryConditions() const {
    return _periodic_boundary;
}

KernelContext::PeriodicBoundaryConditions &KernelContext::periodicBoundaryConditions() {
    return _periodic_boundary;
}

KernelContext::~KernelContext() = default;


}
}






