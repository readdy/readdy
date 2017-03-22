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
 * Topologies are superstructures grouping particles to e.g. molecules.
 *
 * @file Topology.h
 * @brief Definitions for topologies
 * @author clonker
 * @date 26.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <memory>
#include <vector>

#include "BondedPotential.h"
#include "AnglePotential.h"
#include "TorsionPotential.h"
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
    using particles_t = std::vector<std::size_t>;

    Topology(particles_t &&particles);

    Topology(const particles_t &particles);

    Topology(const Topology &) = delete;

    Topology &operator=(const Topology &) = delete;

    Topology(Topology &&) = delete;

    Topology &operator=(Topology &&) = delete;

    virtual ~Topology();

    particles_t::size_type getNParticles() const;

    const particles_t &getParticles() const;

    particles_t &getParticles();

    const std::vector<std::unique_ptr<BondedPotential>> &getBondedPotentials() const;

    const std::vector<std::unique_ptr<AnglePotential>> &getAnglePotentials() const;

    const std::vector<std::unique_ptr<TorsionPotential>> &getTorsionPotentials() const;

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<BondedPotential, T>::value>::type addBondedPotential(Args &&...args) {
        bondedPotentials.push_back(std::make_unique<T>(this, std::forward<Args>(args)...));
    };

    void addBondedPotential(std::unique_ptr<BondedPotential> &&);

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<AnglePotential, T>::value>::type addAnglePotential(Args &&...args) {
        anglePotentials.push_back(std::make_unique<T>(this, std::forward<Args>(args)...));
    };

    void addAnglePotential(std::unique_ptr<AnglePotential> &&);

    template<typename T, typename... Args>
    typename std::enable_if<std::is_base_of<TorsionPotential, T>::value>::type addTorsionPotential(Args &&...args) {
        torsionPotentials.push_back(std::make_unique<T>(this, std::forward<Args>(args)...));
    };

    void addTorsionPotential(std::unique_ptr<TorsionPotential> &&);

    virtual void permuteIndices(const std::vector<std::size_t> &permutation);

protected:
    particles_t particles;
    std::vector<std::unique_ptr<BondedPotential>> bondedPotentials;
    std::vector<std::unique_ptr<AnglePotential>> anglePotentials;
    std::vector<std::unique_ptr<TorsionPotential>> torsionPotentials;
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
