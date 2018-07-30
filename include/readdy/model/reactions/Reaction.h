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
 * Base class for all possible types of reaction. Currently:
 *  - Conversion: A -> B
 *  - Enzymatic: A + C -> B + C where C is a catalyst
 *  - Fission: A -> B + C
 *  - Fusion: A + B -> C
 *
 * @file Reactions.h
 * @brief Reaction base class.
 * @author clonker
 * @date 17.06.16
 */

#pragma once
#include <string>
#include <ostream>
#include <utility>
#include <memory>
#include <spdlog/fmt/ostr.h>
#include <readdy/model/Particle.h>
#include <readdy/model/RandomProvider.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

enum class ReactionType { Conversion, Fusion, Fission, Enzymatic, Decay };

inline std::ostream& operator<<(std::ostream& os, const ReactionType& reactionType) {
    switch (reactionType) {
        case ReactionType::Decay: os << "Decay"; break;
        case ReactionType::Conversion: os << "Conversion"; break;
        case ReactionType::Fusion: os << "Fusion"; break;
        case ReactionType::Fission: os << "Fission"; break;
        case ReactionType::Enzymatic: os << "Enzymatic"; break;
    }
    return os;
}

class Reaction {
public:
    using ReactionId = unsigned short;

    Reaction(std::string name, scalar rate, scalar eductDistance, scalar productDistance, std::uint8_t nEducts,
             std::uint8_t nProducts)
            : _name(std::move(name)), _id(counter++), _rate(rate), _eductDistance(eductDistance),
              _eductDistanceSquared(eductDistance * eductDistance), _productDistance(productDistance),
              _nEducts(nEducts), _nProducts(nProducts) {};

    virtual ~Reaction() = default;

    virtual const ReactionType type() const = 0;

    const std::string &name() const {
        return _name;
    }

    const ReactionId id() const {
        return _id;
    }

    const scalar rate() const {
        return _rate;
    }

    scalar &rate() {
        return _rate;
    }

    std::uint8_t nEducts() const {
        return _nEducts;
    }

    std::uint8_t nProducts() const {
        return _nProducts;
    }

    const scalar eductDistance() const {
        return _eductDistance;
    }

    const scalar eductDistanceSquared() const {
        return _eductDistanceSquared;
    }

    const scalar productDistance() const {
        return _productDistance;
    }

    friend std::ostream &operator<<(std::ostream &os, const Reaction &reaction) {
        os << "Reaction(\"" << reaction._name << "\", N_Educts=" << std::to_string(reaction._nEducts) << ", N_Products="
           << std::to_string(reaction._nProducts) << ", (";
        for (unsigned int i = 0; i < reaction._nEducts; i++) {
            if (i > 0) os << ",";
            os << reaction._educts[i];
        }
        os << ") -> (";
        for (unsigned int i = 0; i < reaction._nProducts; i++) {
            if (i > 0) os << ",";
            os << reaction._products[i];
        }
        os << "), rate=" << reaction._rate << ", eductDist=" << reaction._eductDistance << ", prodDist="
           << reaction._productDistance << ")";
        return os;
    }

    const std::array<ParticleTypeId, 2> &educts() const {
        return _educts;
    };

    const std::array<ParticleTypeId, 2> &products() const {
        return _products;
    };

    const scalar weight1() const {
        return _weight1;
    }

    const scalar weight2() const {
        return _weight2;
    }


protected:
    static ReactionId counter;

    std::uint8_t _nEducts;
    std::uint8_t _nProducts;
    std::array<ParticleTypeId, 2> _educts {{0, 0}};
    std::array<ParticleTypeId, 2> _products {{0, 0}};
    std::string _name;
    ReactionId _id;
    scalar _rate;
    scalar _eductDistance, _eductDistanceSquared;
    scalar _productDistance;

    scalar _weight1 = .5, _weight2 = .5;
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
