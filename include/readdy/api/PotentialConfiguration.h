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
 * This header contains the definitions for bonds, angles, torsion potentials.
 *
 * @file PotentialConfiguration.h
 * @brief Header file containing definitions for topology potential configurations.
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

/**
 * strongly typed enum holding the various bond types
 */
enum class BondType {
    HARMONIC /**< enum value for the classical harmonic spring between two particles */
};
/**
 * strongly typed enum holding the various angle potential types
 */
enum class AngleType {
    HARMONIC /**< enum value for a harmonic potential between a particle triplet fixing a specific angle */
};
/**
 * strongly typed enum holding the various torsion potential types
 */
enum class TorsionType {
    COS_DIHEDRAL /**< enum value for a cosine dihedral potential between a particle quadruple */
};

/**
 * configuration pod class for bonds
 */
struct Bond {
    /**
     * the force constant of the bond, controlling the stiffness
     */
    scalar forceConstant = 0;
    /**
     * the preferred distance w.r.t. this bond
     */
    scalar length = 0;
    /**
     * type of the bond potential
     */
    BondType type = BondType::HARMONIC;
};

/**
 * configuration pod class for angles
 */
struct Angle {
    /**
     * the force constant of the angle, controlling the stiffness
     */
    scalar forceConstant = 0;
    /**
     * the preferred angle between a particle triple w.r.t. this angle potential
     */
    scalar equilibriumAngle = 0;
    /**
     * the type of angle potential
     */
    AngleType type = AngleType::HARMONIC;
};

/**
 * configuration pod class for torsion potentials
 */
struct TorsionAngle {
    /**
     * the stiffness of the torsion potential
     */
    scalar forceConstant {0};
    /**
     * multiplicity of the torsion potential, indicating the number of minima as the bond is rotated through 360 degrees
     */
    scalar multiplicity {0};
    /**
     * the preferred torsion angle w.r.t. this potential
     */
    scalar phi_0 {0};
    /**
     * the type of torsion potential
     */
    TorsionType type = TorsionType::COS_DIHEDRAL;
};

/**
 * potential configuration pod class which serves as a lookup table for particle type tuples/triples/quadruples in order
 * to find their defined potentials
 */
struct PotentialConfiguration {
    /**
     * map that yields a vector of bond configurations for a particle type tuple, map[(a,b)] == map[(b,a)]
     */
    using pair_potential_map = util::particle_type_pair_unordered_map<std::vector<Bond>>;
    /**
     * map that yields a vector of angle configurations for a particle type triple, map[(a,b,c)] == map[(c,b,a)]
     */
    using angle_potential_map = util::particle_type_triple_unordered_map<std::vector<Angle>>;
    /**
     * map that yields a vector of torsion potential configurations, map[(a,b,c,d)] == map[(d,c,b,a)]
     */
    using torsion_potential_map = util::particle_type_quadruple_unordered_map<std::vector<TorsionAngle>>;
    /**
     * a map from particle type tuple to collection of bond potentials
     */
    pair_potential_map pairPotentials;
    /**
     * a map from particle type triple to collection of angle potentials
     */
    angle_potential_map anglePotentials;
    /**
     * a map from particle type quadruple to collection of torsion potentials
     */
    torsion_potential_map torsionPotentials;
};

NAMESPACE_END(api)
NAMESPACE_END(readdy)
