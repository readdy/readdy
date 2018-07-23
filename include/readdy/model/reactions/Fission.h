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
 * This header file contains the declaration of fission reactions, i.e., A->B+C. The products can be placed at a
 * specified distance and weighted by two members.
 *
 * @file Fission.h
 * @brief Declaration of fission reactions, i.e., A->B+C.
 * @author clonker
 * @date 20.06.16
 */

#pragma once

#include <readdy/common/logging.h>
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class Fission : public Reaction {
    using super = Reaction;
public:
    Fission(const std::string &name, ParticleTypeId from, ParticleTypeId to1, ParticleTypeId to2,
            const scalar rate, const scalar productDistance, const scalar weight1 = 0.5,
            const scalar weight2 = 0.5) :
            Reaction(name, rate, 0, productDistance, 1, 2) {
        super::_weight1 = weight1;
        super::_weight2 = weight2;
        _educts = {from};
        _products = {to1, to2};
        const auto sum = weight1 + weight2;
        if (sum != 1) {
            this->_weight1 /= sum;
            this->_weight2 /= sum;
            log::warn("The weights did not add up to 1, they were changed to weight1={}, weight2={}",
                      this->_weight1, this->_weight2);
        }
    }

    const ParticleTypeId getFrom() const {
        return _educts[0];
    }

    const ParticleTypeId getTo1() const {
        return _products[0];
    }

    const ParticleTypeId getTo2() const {
        return _products[1];
    }

    const ReactionType type() const override { return ReactionType::Fission; }

};
NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
