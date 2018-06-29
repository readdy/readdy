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
 * @file Potentials.cpp
 * @brief Core library implementation of order 1 potentials
 * @author clonker
 * @author chrisfroe
 * @date 20.06.16
 */

#include <readdy/model/Kernel.h>
#include <readdy/model/potentials/PotentialsOrder1.h>

namespace readdy {
namespace model {
namespace potentials {

std::string PotentialRegistry::describe() const {
    namespace rus = readdy::util::str;
    auto find_pot_name = [this](ParticleTypeId type) -> const std::string {
        for (auto &&t : _types.get().typeMapping()) {
            if (t.second == type) return t.first;
        }
        return "";
    };
    std::string description;
    if (!potentialsOrder1().empty()) {
        description += fmt::format(" - potentials of order 1:{}", rus::newline);
        for (const auto &types : potentialsOrder1()) {
            description += fmt::format("     * for type \"{}\"{}", find_pot_name(types.first), rus::newline);
            for (auto pot : types.second) {
                description += fmt::format("         * {}{}", pot->describe(), rus::newline);
            }
        }
    }
    if (!potentialsOrder2().empty()) {
        description += fmt::format(" - potentials of order 2:{}", rus::newline);
        for (const auto &types : potentialsOrder2()) {
            description += fmt::format(R"(     * for types "{}" and "{}"{})", find_pot_name(std::get<0>(types.first)),
                                       find_pot_name(std::get<1>(types.first)), rus::newline);
            for (auto pot : types.second) {
                description += fmt::format("         * {}{}", pot->describe(), rus::newline);
            }
        }
    }
    return description;
}

/////////////////////////////////////////////////////////////////////////////
//
// Potentials order 1
//
/////////////////////////////////////////////////////////////////////////////

/**
 * Box potential
 */

Vec3 getMinExtent(const Vec3 &origin, const Vec3 &extent) {
    Vec3 result{0, 0, 0};
    for (auto i = 0; i < 3; i++) {
        if (extent[i] > 0) {
            result[i] = origin[i];
        } else {
            result[i] = origin[i] + extent[i];
        }
    }
    return result;
}

Vec3 getMaxExtent(const Vec3 &origin, const Vec3 &extent) {
    Vec3 result{0, 0, 0};
    for (auto i = 0; i < 3; i++) {
        if (extent[i] > 0) {
            result[i] = origin[i] + extent[i];
        } else {
            result[i] = origin[i];
        }
    }
    return result;
}

Box::Box(particle_type_type particleType, scalar forceConstant, const Vec3 &origin,
                             const Vec3 &extent)
        : super(particleType), origin(origin), extent(extent), forceConstant(forceConstant),
          min(getMinExtent(origin, extent)), max(getMaxExtent(origin, extent)) {}

std::string Box::describe() const {
    return fmt::format("Box potential with origin={}, extent={}, and force constant k={}",
                       origin, extent, forceConstant);
}

std::string Box::type() const {
    return getPotentialName<Box>();
}

/*
 * Sphere Potentials
 */

std::string SphereIn::describe() const {
    return fmt::format("Spherical inclusion potential with origin={}, radius={}, and force constant k={}",
                       origin, radius, forceConstant);
}

std::string SphereIn::type() const {
    return getPotentialName<SphereIn>();
}

std::string SphereOut::describe() const {
    return fmt::format("Spherical exclusion potential with origin={}, radius={}, and force constant k={}",
                       origin, radius, forceConstant);
}

std::string SphereOut::type() const {
    return getPotentialName<SphereOut>();
}

SphericalBarrier::SphericalBarrier(particle_type_type particleType, scalar height, scalar width, const Vec3 &origin, scalar radius)
        : super(particleType), origin(origin), radius(radius), height(height), width(width), r1(radius - width), r2(radius - width / static_cast<scalar>(2.)),
          r3(radius + width / static_cast<scalar>(2.)), r4(radius + width), effectiveForceConstant(static_cast<scalar>(4.) * height / width / width) {
    if (width > radius) {
        throw std::invalid_argument("SphericalBarrier must have a radius larger than its width");
    }
}

std::string SphericalBarrier::describe() const {
    return fmt::format("Spherical barrier potential with origin={}, radius={}, height(energy)={}, and width={}",
                       origin, radius, height, width);
}

std::string SphericalBarrier::type() const {
    return getPotentialName<SphericalBarrier>();
}

/////////////////////////////////////////////////////////////////////////////
//
// Potentials order 2
//
/////////////////////////////////////////////////////////////////////////////

/*
 * Harmonic repulsion
 */

std::string HarmonicRepulsion::describe() const {
    return fmt::format("Harmonic repulsion with force constant k={}", _forceConstant);
}

std::string HarmonicRepulsion::type() const {
    return getPotentialName<HarmonicRepulsion>();
}

/**
 * Weak interaction piecewise harmonic
 */

std::string WeakInteractionPiecewiseHarmonic::describe() const {
    return fmt::format("Weak interaction piecewise harmonic potential with force constant k={}, desired distance={}, "
                               "depth={}, and cutoff={}",
                       forceConstant, conf.desiredParticleDistance, conf.depthAtDesiredDistance,
                       conf.noInteractionDistance);
}

std::string WeakInteractionPiecewiseHarmonic::type() const {
    return getPotentialName<HarmonicRepulsion>();
}

WeakInteractionPiecewiseHarmonic::Configuration::Configuration(const scalar  desiredParticleDistance,
                                                               const scalar  depthAtDesiredDistance,
                                                               const scalar  noInteractionDistance)
        : desiredParticleDistance(desiredParticleDistance), depthAtDesiredDistance(depthAtDesiredDistance),
          noInteractionDistance(noInteractionDistance),
          noInteractionDistanceSquared(noInteractionDistance * noInteractionDistance) {}

LennardJones::LennardJones(particle_type_type type1, particle_type_type type2,
                           unsigned int m, unsigned int n, scalar  cutoffDistance,
                           bool shift, scalar  epsilon, scalar  sigma)
        : super(type1, type2), m(m), n(n),
          cutoffDistance(cutoffDistance), shift(shift), epsilon(epsilon), sigma(sigma),
          cutoffDistanceSquared(cutoffDistance * cutoffDistance) {
    if (m <= n) {
        throw std::invalid_argument("When constructing the LJ potential, the first exponent m=" + std::to_string(m) +
                                    " was not greater than the second exponent n=" + std::to_string(n) + "!");
    }
    auto dm = static_cast<scalar >(m);
    auto dn = static_cast<scalar >(n);
    scalar  r_min = sigma * std::pow(dn / dm, c_::one / (dn - dm));
    k = -epsilon / (std::pow(sigma / r_min, dm) - std::pow(sigma / r_min, dn));
}

std::string LennardJones::describe() const {
    std::string withOrWithout = shift ? "with" : "without";
    return fmt::format("{}-{}-Lennard-Jones potential with cutoff={}, epsilon={}, k={}, and {} energy shift",
                       m, n, cutoffDistance, epsilon, k, withOrWithout);
}

std::string LennardJones::type() const {
    return getPotentialName<LennardJones>();
}

ScreenedElectrostatics::ScreenedElectrostatics(particle_type_type type1, particle_type_type type2,
                                               scalar  electrostaticStrength, scalar  inverseScreeningDepth,
                                               scalar  repulsionStrength, scalar  repulsionDistance, unsigned int exponent,
                                               scalar  cutoff)
        : super(type1, type2), electrostaticStrength(electrostaticStrength), inverseScreeningDepth(inverseScreeningDepth),
          repulsionStrength(repulsionStrength), repulsionDistance(repulsionDistance), exponent(exponent), cutoff(cutoff),
          cutoffSquared(cutoff * cutoff) {
    if (inverseScreeningDepth < 0) {
        throw std::invalid_argument("inverse screening depth must be positive!");
    }
    if (repulsionStrength < 0) {
        throw std::invalid_argument("repulsion strength must be positive!");
    }
    if (repulsionDistance < 0) {
        throw std::invalid_argument("repulsion distance must be positive!");
    }
    if (cutoff < 0) {
        throw std::invalid_argument("cutoff must be positive!");
    }
}

std::string ScreenedElectrostatics::describe() const {
    return fmt::format("Screened electrostatics potential with electrostatic strength={}, inverse screening depth={}, "
                               "repulsion strength={}, repulsion distance={}, exponent={}, and cutoff={}",
                       electrostaticStrength, inverseScreeningDepth, repulsionStrength, repulsionDistance,
                       exponent, cutoff);
}

std::string ScreenedElectrostatics::type() const {
    return getPotentialName<ScreenedElectrostatics>();
}

}
}
}
