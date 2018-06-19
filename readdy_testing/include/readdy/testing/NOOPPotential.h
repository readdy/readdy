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
