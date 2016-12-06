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
 * @file SingleCPUReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 22.06.16
 */

#include <readdy/model/Particle.h>
#include <readdy/kernel/singlecpu/reactions/SCPUReactions.h>

using rdy_particle_t = readdy::model::Particle;

namespace readdy {
namespace kernel {
namespace scpu {
namespace reactions {
void SCPUConversion::perform(const rdy_particle_t &p1_in, const rdy_particle_t &p2_in, rdy_particle_t &p1_out,
                         rdy_particle_t &p2_out, rnd_normal rnd) const {
    p1_out.setPos(p1_in.getPos());
    p1_out.setType(getTypeTo());
    p1_out.setId(p1_in.getId());
}

void SCPUEnzymatic::perform(const rdy_particle_t &p1_in, const rdy_particle_t &p2_in, rdy_particle_t &p1_out,
                        rdy_particle_t &p2_out, rnd_normal rnd) const {
    if (p1_in.getType() == getCatalyst()) {
        // p1 is the catalyst
        p1_out.setType(getCatalyst());
        p1_out.setPos(p1_in.getPos());
        p1_out.setId(p1_in.getId());
        p2_out.setType(getTo());
        p2_out.setPos(p2_in.getPos());
    } else {
        // p2 is the catalyst
        p1_out.setType(getCatalyst());
        p1_out.setPos(p2_in.getPos());
        p1_out.setId(p2_in.getId());
        p2_out.setType(getTo());
        p2_out.setPos(p1_in.getPos());
    }
}

void SCPUFission::perform(const rdy_particle_t &p1_in, const rdy_particle_t &p2_in, rdy_particle_t &p1_out,
                      rdy_particle_t &p2_out, rnd_normal rnd) const {
    // as long as the orientation is uniform, it does not matter of which type p1_in and p2_in are.
    auto n3 = rnd(0, 1);
    n3 /= std::sqrt(n3 * n3);
    p1_out.setType(getTo1());
    p1_out.setPos(p1_in.getPos() + getWeight1() * getProductDistance() * n3);

    p2_out.setType(getTo2());
    p2_out.setPos(p1_in.getPos() - getWeight2() * getProductDistance() * n3);
}

void SCPUFusion::perform(const rdy_particle_t &p1_in, const rdy_particle_t &p2_in, rdy_particle_t &p1_out,
                     rdy_particle_t &p2_out, rnd_normal rnd) const {
    p1_out.setType(getTo());
    if (getFrom1() == p1_in.getType()) {
        p1_out.setPos(p1_in.getPos() + getWeight1() * (p2_in.getPos() - p1_in.getPos()));
    } else {
        p1_out.setPos(p2_in.getPos() + getWeight1() * (p1_in.getPos() - p2_in.getPos()));
    }
}

}
}
}
}



