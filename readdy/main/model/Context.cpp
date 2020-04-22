/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * @file KernelContext.cpp
 * @brief Implementation file of the KernelContext.
 * @author clonker
 * @author chrisfroe
 * @date 18.04.16
 * @todo make proper reference to KernelContext.h, is kBT really indepdendent of t?
 */

#include <readdy/model/Context.h>

#include <readdy/api/KernelConfiguration.h>

#include <readdy/common/Utils.h>
#include <readdy/model/_internal/Util.h>
#include <readdy/common/boundary_condition_operations.h>
#include <readdy/model/potentials/PotentialsOrder1.h>

namespace readdy::model {

using particle_t = readdy::model::Particle;

void Context::setKernelConfiguration(const std::string &s) {
    _kernelConfiguration = nlohmann::json::parse(s);
}

namespace {
bool boxPotentialValid(const Context &self, potentials::Box *potential) {
    auto bbox = self.getBoxBoundingVertices();
    auto pbc = self.periodicBoundaryConditions();
    bool valid = true;
    for(std::size_t dim = 0; dim < 3; ++dim) {
        auto periodic = pbc.at(dim);
        auto bbox0Pos = std::get<0>(bbox)[dim];
        auto bbox1Pos = std::get<1>(bbox)[dim];
        auto potential1Pos = potential->getOrigin()[dim];
        auto potential2Pos = potential1Pos + potential->getExtent()[dim];
        if(!periodic) {
            valid &= bbox0Pos < potential1Pos && bbox1Pos > potential2Pos;
        }
    }
    return valid;
}
}

void Context::validate() const {
    auto periodic = std::accumulate(periodicBoundaryConditions().begin(),
                                    periodicBoundaryConditions().end(), true, std::logical_and<>());
    if(!periodic) {
        // check if there are box potentials for each particle type and that these box potentials are valid
        for(const auto &entry : particleTypes().typeMapping()) {
            auto ptype = entry.second;
            auto potIt = potentials().potentialsOrder1().find(ptype);
            bool valid = true;
            if(potIt != potentials().potentialsOrder1().end()) {
                bool gotValidBoxPotential = false;
                for(const auto potPtr : potIt->second) {
                    if(potPtr->type() == potentials::getPotentialName<potentials::Box>()) {
                        gotValidBoxPotential |= boxPotentialValid(*this, dynamic_cast<potentials::Box*>(potPtr));
                    }
                }
                valid &= gotValidBoxPotential;
            } else {
                valid = false;
            }
            if(!valid) {
                throw std::logic_error(fmt::format("For particle type {} there was no valid box potential in direction "
                                                           "of the non-periodic boundaries configured.", entry.first));
            }
        }
    }
}

std::string Context::describe() {
    std::string description;
    description += fmt::format("Configured kernel context with:\n");
    description += fmt::format("--------------------------------\n");
    description += fmt::format(" - kBT = {}\n", kBT());
    description += fmt::format(" - periodic b.c. = ({}, {}, {})\n", periodicBoundaryConditions()[0],
            periodicBoundaryConditions()[1], periodicBoundaryConditions()[2]);
    description += fmt::format(" - box size = ({}, {}, {})\n", boxSize()[0], boxSize()[1], boxSize()[2]);

    description += _particleTypeRegistry.describe();
    description += _potentialRegistry.describe();
    description += _reactionRegistry.describe();
    description += _topologyRegistry.describe();
    return description;
}

scalar Context::calculateMaxCutoff() const {
    scalar max_cutoff{0};
    for (const auto &entry : potentials().potentialsOrder2()) {
        for (const auto &potential : entry.second) {
            max_cutoff = std::max(max_cutoff, potential->getCutoffRadius());
        }
    }
    for (const auto &entry : reactions().order2()) {
        for (const auto &reaction : entry.second) {
            max_cutoff = std::max(max_cutoff, reaction->eductDistance());
        }
    }
    for (const auto &entry : _topologyRegistry.spatialReactionRegistry()) {
        for (const auto &reaction : entry.second) {
            max_cutoff = std::max(max_cutoff, reaction.radius());
        }
    }
    return max_cutoff;
}

void Context::setTypeRegistryReferences() {
    _reactionRegistry._types = &_particleTypeRegistry;
    _potentialRegistry._types = &_particleTypeRegistry;
    _compartmentRegistry._types = &_particleTypeRegistry;
    _topologyRegistry._types = &_particleTypeRegistry;
}

Context::Context() {
    setTypeRegistryReferences();
}

Context::Context(const Context &rhs) { *this = rhs; }

Context::Context(Context &&rhs) noexcept { *this = std::move(rhs); }

Context &Context::operator=(Context &&rhs) noexcept {
    _particleTypeRegistry = std::move(rhs._particleTypeRegistry);
    _reactionRegistry = std::move(rhs._reactionRegistry);
    _potentialRegistry = std::move(rhs._potentialRegistry);
    _compartmentRegistry = std::move(rhs._compartmentRegistry);
    _topologyRegistry = std::move(rhs._topologyRegistry);
    _kernelConfiguration = rhs._kernelConfiguration;
    _kBT = rhs._kBT;
    _box_size = rhs._box_size;
    _periodic_boundary = rhs._periodic_boundary;
    _recordReactionsWithPositions = rhs._recordReactionsWithPositions;
    _recordReactionCounts = rhs._recordReactionCounts;
    _recordVirial = rhs._recordVirial;
    setTypeRegistryReferences();
    return *this;
}

Context &Context::operator=(const Context &rhs) {
    if (this != &rhs) {
        _particleTypeRegistry = rhs._particleTypeRegistry;
        _reactionRegistry = rhs._reactionRegistry;
        _potentialRegistry = rhs._potentialRegistry;
        _compartmentRegistry = rhs._compartmentRegistry;
        _topologyRegistry = rhs._topologyRegistry;
        _kernelConfiguration = rhs._kernelConfiguration;
        _kBT = rhs._kBT;
        _box_size = rhs._box_size;
        _periodic_boundary = rhs._periodic_boundary;
        _recordReactionsWithPositions = rhs._recordReactionsWithPositions;
        _recordReactionCounts = rhs._recordReactionCounts;
        _recordVirial = rhs._recordVirial;
        setTypeRegistryReferences();
    }
    return *this;
}

}
