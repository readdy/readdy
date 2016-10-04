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

using _rdy_particle_t = readdy::model::Particle;

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace reactions {
void Conversion::perform(const _rdy_particle_t &p1_in, const _rdy_particle_t &p2_in, _rdy_particle_t &p1_out,
                         _rdy_particle_t &p2_out, const rnd_ptr &rnd) const {
    p1_out.setPos(p1_in.getPos());
    p1_out.setType(getTypeTo());
    p1_out.setId(p1_in.getId());
}

Conversion *Conversion::replicate() const {
    return new Conversion(*this);
}


void Enzymatic::perform(const _rdy_particle_t &p1_in, const _rdy_particle_t &p2_in, _rdy_particle_t &p1_out,
                        _rdy_particle_t &p2_out, const rnd_ptr &rnd) const {
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

Enzymatic *Enzymatic::replicate() const {
    return new Enzymatic(*this);
}

void Fission::perform(const _rdy_particle_t &p1_in, const _rdy_particle_t &p2_in, _rdy_particle_t &p1_out,
                      _rdy_particle_t &p2_out, const rnd_ptr &rnd) const {
    // as long as the orientation is uniform, it does not matter of which type p1_in and p2_in are.
    auto n3 = rnd->getNormal3();
    n3 /= sqrt(n3 * n3);
    p1_out.setType(getTo1());
    p1_out.setPos(p1_in.getPos() + getWeight1() * getProductDistance() * n3);

    p2_out.setType(getTo2());
    p2_out.setPos(p1_in.getPos() - getWeight2() * getProductDistance() * n3);
}

Fission *Fission::replicate() const {
    return new Fission(*this);
}

void Fusion::perform(const _rdy_particle_t &p1_in, const _rdy_particle_t &p2_in, _rdy_particle_t &p1_out,
                     _rdy_particle_t &p2_out, const rnd_ptr &rnd) const {
    p1_out.setType(getTo());
    if (getFrom1() == p1_in.getType()) {
        p1_out.setPos(p1_in.getPos() + getWeight1() * (p2_in.getPos() - p1_in.getPos()));
    } else {
        p1_out.setPos(p2_in.getPos() + getWeight1() * (p1_in.getPos() - p2_in.getPos()));
    }
}

Fusion *Fusion::replicate() const {
    return new Fusion(*this);
}
}
}
}
}



