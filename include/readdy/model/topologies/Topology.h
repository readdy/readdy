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

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

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
    using particle_indices = std::vector<std::size_t>;
    using particle_index = particle_indices::value_type;

    using bonded_potential = pot::BondedPotential;
    using harmonic_bond = TopologyActionFactory::harmonic_bond;

    using angle_potential = pot::AnglePotential;
    using harmonic_angle = TopologyActionFactory::harmonic_angle;

    using torsion_potential = pot::TorsionPotential;
    using cos_dihedral = TopologyActionFactory::cos_dihedral;

    explicit Topology(particle_indices particles) : particles(std::move(particles)) { }

    Topology(const Topology &) = delete;

    Topology &operator=(const Topology &) = delete;

    Topology(Topology &&) = default;

    Topology &operator=(Topology &&) = default;

    virtual ~Topology() = default;

    particle_indices::size_type getNParticles() const {
        return particles.size();
    }

    const particle_indices &getParticles() const {
        return particles;
    }

    particle_indices &getParticles() {
        return particles;
    }

    const std::vector<std::unique_ptr<bonded_potential>> &getBondedPotentials() const {
        return bondedPotentials;
    }

    const std::vector<std::unique_ptr<angle_potential>> &getAnglePotentials() const {
        return anglePotentials;
    }

    const std::vector<std::unique_ptr<torsion_potential>> &getTorsionPotentials() const {
        return torsionPotentials;
    }

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<bonded_potential, T>::value>::type addBondedPotential(Args &&...args) {
        bondedPotentials.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    }

    void addBondedPotential(std::unique_ptr<bonded_potential> &&pot) {
        const auto n = getNParticles();
        for(const auto& bond : pot->getBonds()) {
            if (bond.idx1 >= n) {
                throw std::invalid_argument("the first particle (" + std::to_string(bond.idx1) + ") was out of bounds!");
            }
            if (bond.idx2 >= n) {
                throw std::invalid_argument("the second particle (" + std::to_string(bond.idx2) + ") was out of bounds!");
            }
        }
        bondedPotentials.push_back(std::move(pot));
    }

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<angle_potential, T>::value>::type addAnglePotential(Args &&...args) {
        anglePotentials.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    }

    void addAnglePotential(std::unique_ptr<angle_potential> &&pot) {
        anglePotentials.push_back(std::move(pot));
    }

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<torsion_potential, T>::value>::type addTorsionPotential(Args &&...args) {
        torsionPotentials.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    }

    void addTorsionPotential(std::unique_ptr<torsion_potential> &&pot) {
        torsionPotentials.push_back(std::move(pot));
    }

    virtual void permuteIndices(const std::vector<std::size_t> &permutation) {
        std::transform(particles.begin(), particles.end(), particles.begin(), [&permutation](std::size_t index) {
            return permutation[index];
        });
    }

protected:
    particle_indices particles;
    std::vector<std::unique_ptr<bonded_potential>> bondedPotentials;
    std::vector<std::unique_ptr<angle_potential>> anglePotentials;
    std::vector<std::unique_ptr<torsion_potential>> torsionPotentials;
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
