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
 * @file BondPotential.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/BondedPotential.h>
#include <readdy/model/topologies/Topology.h>

namespace readdy {
namespace model {
namespace top {

/*
 * Super class
 */

BondedPotential::BondedPotential(Topology *const topology) : TopologyPotential(topology) {}

/*
 * Harmonic bond
 */

HarmonicBondPotential::HarmonicBondPotential(Topology *const topology, const std::vector<Bond> &bonds)
        : BondedPotential(topology), bonds(bonds) {}
HarmonicBondPotential::HarmonicBondPotential(Topology *const topology, std::vector<Bond> bonds)
        : BondedPotential(topology), bonds(std::move(bonds)) {}

void HarmonicBondPotential::addBond(std::size_t idx1, std::size_t idx2, double length, double forceConstant) {
    const auto n = topology->getNParticles();
    if(idx1 >= n) {
        throw std::invalid_argument("the first particle index (" + std::to_string(idx1) + ") was out of bounds!");
    }
    if(idx2 >= n) {
        throw std::invalid_argument("the second particle index (" + std::to_string(idx2) + ") was out of bounds!");
    }
    bonds.emplace_back(idx1, idx2, length, forceConstant);
}


}
}
}