/**
 * << detailed description >>
 *
 * @file SingleCPUReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 22.06.16
 */

#include <readdy/model/Particle.h>
#include <readdy/kernel/singlecpu/reactions/SingleCPUReactions.h>

using rdy_particle_t = readdy::model::Particle;

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace reactions {
void Conversion::perform(const rdy_particle_t &p1_in, const rdy_particle_t &p2_in, rdy_particle_t &p1_out,
                         rdy_particle_t &p2_out, rnd_normal rnd) const {
    p1_out.setPos(p1_in.getPos());
    p1_out.setType(getTypeTo());
    p1_out.setId(p1_in.getId());
}

void Enzymatic::perform(const rdy_particle_t &p1_in, const rdy_particle_t &p2_in, rdy_particle_t &p1_out,
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

void Fission::perform(const rdy_particle_t &p1_in, const rdy_particle_t &p2_in, rdy_particle_t &p1_out,
                      rdy_particle_t &p2_out, rnd_normal rnd) const {
    // as long as the orientation is uniform, it does not matter of which type p1_in and p2_in are.
    auto n3 = rnd(0, 1);
    n3 /= std::sqrt(n3 * n3);
    p1_out.setType(getTo1());
    p1_out.setPos(p1_in.getPos() + getWeight1() * getProductDistance() * n3);

    p2_out.setType(getTo2());
    p2_out.setPos(p1_in.getPos() - getWeight2() * getProductDistance() * n3);
}

void Fusion::perform(const rdy_particle_t &p1_in, const rdy_particle_t &p2_in, rdy_particle_t &p1_out,
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



