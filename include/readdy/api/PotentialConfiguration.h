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
NAMESPACE_BEGIN(api)

enum class BondType { HARMONIC };
enum class AngleType { HARMONIC };
enum class TorsionType { COS_DIHEDRAL };

struct Bond {
    scalar forceConstant = 0, length = 0;
    BondType type = BondType::HARMONIC;
};

struct Angle {
    scalar forceConstant = 0, equilibriumAngle = 0;
    AngleType type = AngleType::HARMONIC;
};

struct TorsionAngle {
    scalar forceConstant {0}, multiplicity {0}, phi_0 {0};
    TorsionType type = TorsionType::COS_DIHEDRAL;
};

struct PotentialConfiguration {
    using pair_potential_map = util::particle_type_pair_unordered_map<std::vector<Bond>>;
    using angle_potential_map = util::particle_type_triple_unordered_map<std::vector<Angle>>;
    using torsion_potential_map = util::particle_type_quadruple_unordered_map<std::vector<TorsionAngle>>;
    pair_potential_map pairPotentials;
    angle_potential_map anglePotentials;
    torsion_potential_map torsionPotentials;
};

NAMESPACE_END(api)
NAMESPACE_END(readdy)
