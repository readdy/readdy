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
 * Second order potential which is only useful for testing. The behavior is defined by the parameters cutoff, force and energy
 * which the constructor accepts. This approach is similar to mocking, but actually mocking a PotentialOrder2 has too many side effects.
 *
 * @file NOOPPotential.h
 * @brief Second-order potential that has no real effect. Just returns values defined during construction.
 * @author chrisfroe
 * @date 19.08.16
 */

#ifndef READDY_MAIN_NOOPPOTENTIAL_H
#define READDY_MAIN_NOOPPOTENTIAL_H

#include <readdy/model/potentials/PotentialOrder2.h>

namespace readdy {
namespace testing {
struct NOOPPotentialOrder2 : public readdy::model::potentials::PotentialOrder2 {
    NOOPPotentialOrder2(particle_type_type particleType1, particle_type_type particleType2,
                        readdy::scalar cutoff = 0, readdy::scalar force = 0, readdy::scalar energy = 0)
            : PotentialOrder2(particleType1, particleType2), cutoff(cutoff), force(force), energy(energy) {}

    std::string describe() const override {
        return "NOOPPotential with types " + std::to_string(_particleType1) + ", " + std::to_string(_particleType2);
    }

    readdy::scalar calculateEnergy(const Vec3 &x_ij) const override {
        return energy;
    }

    std::string type() const override {
        return "noop pot";
    }

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override {
        force[0] = NOOPPotentialOrder2::force;
        force[1] = NOOPPotentialOrder2::force;
        force[2] = NOOPPotentialOrder2::force;
    }

    readdy::scalar getCutoffRadiusSquared() const override {
        return cutoff*cutoff;
    }

    readdy::scalar cutoff, force, energy;
};
}
}

#endif //READDY_MAIN_NOOPPOTENTIAL_H
