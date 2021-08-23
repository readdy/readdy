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
 * @file Potentials.cpp
 * @brief Core library implementation of order 1 potentials
 * @author clonker
 * @author chrisfroe
 * @date 20.06.16
 */

#include <readdy/model/Kernel.h>
#include <readdy/model/potentials/PotentialsOrder1.h>

namespace readdy::model::potentials {

std::string PotentialRegistry::describe() const {
    namespace rus = readdy::util::str;
    auto find_pot_name = [this](ParticleTypeId type) -> std::string {
        for (auto &&t : _types->typeMapping()) {
            if (t.second == type) return t.first;
        }
        return "";
    };
    std::string description;
    if (!potentialsOrder1().empty()) {
        description += fmt::format(" - potentials of order 1:\n");
        for (const auto &types : potentialsOrder1()) {
            description += fmt::format("     * for type \"{}\"\n", find_pot_name(types.first));
            for (auto pot : types.second) {
                description += fmt::format("         * {}\n", pot->describe());
            }
        }
    }
    if (!potentialsOrder2().empty()) {
        description += fmt::format(" - potentials of order 2:\n");
        for (const auto &types : potentialsOrder2()) {
            description += fmt::format(R"(     * for types "{}" and "{}"{})", find_pot_name(std::get<0>(types.first)),
                                       find_pot_name(std::get<1>(types.first)), "\n");
            for (auto pot : types.second) {
                description += fmt::format("         * {}\n", pot->describe());
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

SphericalBarrier::SphericalBarrier(ParticleTypeId particleType, scalar height, scalar width, const Vec3 &origin, scalar radius)
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

/**
 * Cylindrical potentials
 */

template<>
std::string Cylinder<true>::type() const {
    return getPotentialName<Cylinder<true>>();
}

template<>
std::string Cylinder<false>::type() const {
    return getPotentialName<Cylinder<false>>();
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
    return fmt::format("Harmonic repulsion with Force constant k={}", _forceConstant);
}

std::string HarmonicRepulsion::type() const {
    return getPotentialName<HarmonicRepulsion>();
}

/**
 * Weak interaction piecewise harmonic
 */

std::string WeakInteractionPiecewiseHarmonic::describe() const {
    return fmt::format("Weak interaction piecewise harmonic potential with Force constant k={}, desired distance={}, "
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

LennardJones::LennardJones(ParticleTypeId type1, ParticleTypeId type2,
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
    scalar  r_min = sigma * std::pow(dn / dm, 1. / (dn - dm));
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

ScreenedElectrostatics::ScreenedElectrostatics(ParticleTypeId type1, ParticleTypeId type2,
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
