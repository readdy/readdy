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
 * @file PotentialConfiguration.h
 * @brief << brief description >>
 * @author clonker
 * @date 17.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <unordered_map>

#include <readdy/common/common.h>
#include <readdy/common/ParticleTypeTuple.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(graph)

enum class BondType { HARMONIC };
enum class AngleType { HARMONIC };
enum class TorsionType { COS_DIHEDRAL };

struct Bond {
    double forceConstant, length;
    BondType type = BondType::HARMONIC;
};

struct Angle {
    double forceConstant, equilibriumAngle;
    AngleType type = AngleType::HARMONIC;
};

struct TorsionAngle {
    double forceConstant, phi_0, multiplicity;
    TorsionType type = TorsionType::COS_DIHEDRAL;
};

struct PotentialConfiguration {
    using pair_potential_map = std::unordered_map<util::particle_type_pair, std::vector<Bond>, util::particle_type_pair_hasher, util::particle_type_pair_equal_to>;
    using angle_potential_map = std::unordered_map<util::particle_type_triple, std::vector<Angle>, util::particle_type_triple_hasher, util::particle_type_triple_equal_to>;
    using torsion_potential_map = std::unordered_map<util::particle_type_quadruple, std::vector<TorsionAngle>, util::particle_type_quadruple_hasher, util::particle_type_quadruple_equal_to>;
    pair_potential_map pairPotentials;
    angle_potential_map anglePotentials;
    torsion_potential_map torsionPotentials;
};

NAMESPACE_END(graph)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
