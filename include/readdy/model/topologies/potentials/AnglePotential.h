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
 * << detailed description >>
 *
 * @file AnglePotential.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.01.17
 * @copyright BSD-3
 */

#pragma once

#include <cstddef>
#include <tuple>
#include <vector>
#include "TopologyPotential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(pot)

class AnglePotential : public TopologyPotential {
public:
    using angles = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;

    AnglePotential() : TopologyPotential() {}

    ~AnglePotential() override = default;
};

struct AngleConfiguration {

    AngleConfiguration(size_t idx1, size_t idx2, size_t idx3, scalar forceConstant, scalar theta_0)
            : idx1(idx1), idx2(idx2), idx3(idx3), equilibriumAngle(theta_0), forceConstant(forceConstant) {}

    const std::size_t idx1, idx2, idx3;
    const scalar equilibriumAngle, forceConstant;
};


class HarmonicAnglePotential : public AnglePotential {
public:
    using angle = AngleConfiguration;
    using angle_configurations = std::vector<AngleConfiguration>;

    explicit HarmonicAnglePotential(angle_configurations angles) : AnglePotential(), angles(std::move(angles)) {}
    HarmonicAnglePotential(const HarmonicAnglePotential&) = default;
    HarmonicAnglePotential& operator=(const HarmonicAnglePotential&) = delete;
    HarmonicAnglePotential(HarmonicAnglePotential&&) = default;
    HarmonicAnglePotential& operator=(HarmonicAnglePotential&&) = delete;

    ~HarmonicAnglePotential() override = default;

    std::unique_ptr<EvaluatePotentialAction>
    createForceAndEnergyAction(const TopologyActionFactory *factory) override;

    const angle_configurations &getAngles() const {
        return angles;
    }

    scalar calculateEnergy(const Vec3 &x_ji, const Vec3 &x_jk, const angle &angle) const;

    void
    calculateForce(Vec3 &f_i, Vec3 &f_j, Vec3 &f_k, const Vec3 &x_ji, const Vec3 &x_jk, const angle &angle) const;

protected:
    angle_configurations angles;
};

NAMESPACE_END(pot)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
