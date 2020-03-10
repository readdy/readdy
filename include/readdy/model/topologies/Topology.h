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
 * Topologies are superstructures grouping particles to e.g. molecules.
 *
 * @file Topology.h
 * @brief Definitions for topologies
 * @author clonker
 * @date 26.01.17
 * @copyright BSD-3
 */

#pragma once

#include <memory>
#include <vector>

#include "TopologyActionFactory.h"

namespace readdy::model::top {

/**
 * Topologies are superstructures grouping particles to e.g. molecules.
 *
 * A topology consists of:
 *  - particle indices
 *  - potentials between particles (bonds, angles, torsions)
 *
 *  It is created by specifying a set of particles and creating corresponding potentials between these.
 */
class Topology {
public:
    using BondedPotential = pot::BondedPotential;
    using HarmonicBond = TopologyActionFactory::harmonic_bond;

    using AnglePotential = pot::AnglePotential;
    using HarmonicAngle = TopologyActionFactory::harmonic_angle;

    using TorsionPotential = pot::TorsionPotential;
    using CosineDihedral = TopologyActionFactory::cos_dihedral;

    Topology() = default;

    Topology(const Topology &) = delete;

    Topology &operator=(const Topology &) = delete;

    Topology(Topology &&) = default;

    Topology &operator=(Topology &&) = default;

    virtual ~Topology() = default;

    [[nodiscard]] const std::vector<std::unique_ptr<BondedPotential>> &getBondedPotentials() const {
        return bondedPotentials;
    }

    [[nodiscard]] const std::vector<std::unique_ptr<AnglePotential>> &getAnglePotentials() const {
        return anglePotentials;
    }

    [[nodiscard]] const std::vector<std::unique_ptr<TorsionPotential>> &getTorsionPotentials() const {
        return torsionPotentials;
    }

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<BondedPotential, T>::value>::type addBondedPotential(Args &&...args) {
        bondedPotentials.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    }

    void addBondedPotential(std::unique_ptr<BondedPotential> &&pot) {
        bondedPotentials.push_back(std::move(pot));
    }

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<AnglePotential, T>::value>::type addAnglePotential(Args &&...args) {
        anglePotentials.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    }

    void addAnglePotential(std::unique_ptr<AnglePotential> &&pot) {
        anglePotentials.push_back(std::move(pot));
    }

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<TorsionPotential, T>::value>::type addTorsionPotential(Args &&...args) {
        torsionPotentials.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    }

    void addTorsionPotential(std::unique_ptr<TorsionPotential> &&pot) {
        torsionPotentials.push_back(std::move(pot));
    }

protected:
    std::vector<std::unique_ptr<BondedPotential>> bondedPotentials;
    std::vector<std::unique_ptr<AnglePotential>> anglePotentials;
    std::vector<std::unique_ptr<TorsionPotential>> torsionPotentials;
};

}
