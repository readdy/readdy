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
 * @file Bond.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.01.17
 * @copyright GPL-3
 */

#pragma once
#include <cstddef>
#include <tuple>
#include <vector>
#include <string>
#include "TopologyPotential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(pot)

struct BondConfiguration {
    BondConfiguration(std::size_t idx1, std::size_t idx2, scalar forceConstant, scalar length)
            : idx1(idx1), idx2(idx2), length(length), forceConstant(forceConstant) {}

    std::size_t idx1, idx2;
    scalar length, forceConstant;
};


class BondedPotential : public TopologyPotential {
public:
    using bond_configuration = BondConfiguration;
    using bond_configurations = std::vector<bond_configuration>;

    explicit BondedPotential(bond_configurations bonds) : TopologyPotential(), bonds(std::move(bonds)) {}

    BondedPotential(const BondedPotential&) = default;
    BondedPotential& operator=(const BondedPotential&) = delete;
    BondedPotential(BondedPotential&&) = default;
    BondedPotential& operator=(BondedPotential&&) = delete;
    ~BondedPotential() override = default;

    const bond_configurations &getBonds() const {
        return bonds;
    }
protected:
    bond_configurations bonds;
};


class HarmonicBondPotential : public BondedPotential {
public:

    explicit HarmonicBondPotential(const bond_configurations &bonds) : BondedPotential(bonds) {};
    HarmonicBondPotential(const HarmonicBondPotential&) = default;
    HarmonicBondPotential& operator=(const HarmonicBondPotential&) = delete;
    HarmonicBondPotential(HarmonicBondPotential&&) = default;
    HarmonicBondPotential& operator=(HarmonicBondPotential&&) = delete;

    ~HarmonicBondPotential() override = default;

    scalar calculateEnergy(const Vec3 &x_ij, const bond_configuration &bond) const {
        const auto norm = std::sqrt(x_ij * x_ij);
        return bond.forceConstant * (norm - bond.length) * (norm - bond.length);
    }

    void calculateForce(Vec3 &force, const Vec3 &x_ij, const bond_configuration &bond) const {
        const auto norm = x_ij.norm();
        force += (2. * bond.forceConstant * (norm - bond.length) / norm) * x_ij;
    }

    std::unique_ptr<EvaluatePotentialAction> createForceAndEnergyAction(const TopologyActionFactory *) override;

};

NAMESPACE_END(pot)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
