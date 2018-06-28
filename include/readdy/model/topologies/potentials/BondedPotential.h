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
 * @file Bond.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.01.17
 * @copyright GNU Lesser General Public License v3.0
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
