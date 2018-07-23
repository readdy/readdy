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
 * Declaration of fusion reactions, i.e., A+B->C. Additionally to the types, they also have a rate, an educt
 * distance and two weights, which determine where (on a line between A and B) C should be placed.
 *
 * @file Fusion.h
 * @brief Header file containing the declaration for fusion reactions, i.e., A+B->C.
 * @author clonker
 * @date 20.06.16
 */

#pragma once
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class Fusion : public Reaction {
    using super = Reaction;
public:
    Fusion(const std::string &name, ParticleTypeId from1, ParticleTypeId from2, ParticleTypeId to,
           const scalar rate, const scalar eductDistance, const scalar weight1 = 0.5,
           const scalar weight2 = 0.5) : Reaction(name, rate, eductDistance, 0, 2, 1){
        super::_weight1 = weight1;
        super::_weight2 = weight2;
        _educts = {from1, from2};
        _products = {to};

        const auto sum = weight1 + weight2;
        if (sum != 1) {
            this->_weight1 /= sum;
            this->_weight2 /= sum;
            log::warn("The weights did not add up to 1, they were changed to weight1={}, weight2={}",
                                 this->_weight1, this->_weight2);
        }
    }

    const ParticleTypeId getFrom1() const {
        return _educts[0];
    }

    const ParticleTypeId getFrom2() const {
        return _educts[1];
    }

    const ParticleTypeId getTo() const {
        return _products[0];
    }

    const ReactionType type() const override {
        return ReactionType::Fusion;
    }


};
NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
