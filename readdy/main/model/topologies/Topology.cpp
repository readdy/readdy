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
 * @file Topology.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/Topology.h>

namespace readdy {
namespace model {
namespace top {
readdy::model::top::Topology::~Topology() = default;

Topology::Topology(Topology::particles_t &&p) : particles(std::move(p)) {}

Topology::particles_t::size_type Topology::getNParticles() const {
    return particles.size();
}

const Topology::particles_t &Topology::getParticles() const {
    return particles;
}

const std::vector<std::unique_ptr<BondedPotential>> &Topology::getBondedPotentials() const {
    return bondedPotentials;
}

const std::vector<std::unique_ptr<AnglePotential>> &Topology::getAnglePotentials() const {
    return anglePotentials;
}

const std::vector<std::unique_ptr<TorsionPotential>> &Topology::getTorsionPotentials() const {
    return torsionPotentials;
}

void Topology::addBondedPotential(std::unique_ptr<BondedPotential> &&pot) {
    bondedPotentials.push_back(std::move(pot));
}

void Topology::addAnglePotential(std::unique_ptr<AnglePotential> &&pot) {
    anglePotentials.push_back(std::move(pot));
}

void Topology::addTorsionPotential(std::unique_ptr<TorsionPotential> &&pot) {
    torsionPotentials.push_back(std::move(pot));
}

}
}
}
