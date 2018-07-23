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
 * @file TopologyPotentialActionFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.04.17
 * @copyright GPL-3
 */

#pragma once

#include <readdy/common/macros.h>
#include "TopologyPotentialActions.h"
#include "BondedPotential.h"
#include "AnglePotential.h"
#include "TorsionPotential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(pot)

class TopologyPotentialActionFactory {
public:

    using harmonic_bond = CalculateHarmonicBondPotential::harmonic_bond;
    using harmonic_angle = CalculateHarmonicAnglePotential::harmonic_angle;
    using cos_dihedral = CalculateCosineDihedralPotential::cos_dihedral;

    virtual std::unique_ptr<CalculateHarmonicBondPotential>
    createCalculateHarmonicBondPotential(const harmonic_bond *) const = 0;

    virtual std::unique_ptr<CalculateHarmonicAnglePotential>
    createCalculateHarmonicAnglePotential(const harmonic_angle *) const = 0;

    virtual std::unique_ptr<CalculateCosineDihedralPotential>
    createCalculateCosineDihedralPotential(const cos_dihedral *) const = 0;
};

NAMESPACE_END(pot)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
