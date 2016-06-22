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
#include <math.h>

using particle_t = readdy::model::Particle;

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace reactions {
                void SingleCPUConversion::perform(const particle_t &p1_in, const particle_t &p2_in, particle_t &p1_out, particle_t &p2_out) const {
                    p1_out.setPos(p1_in.getPos());
                    p1_out.setType(getTypeTo());
                    p1_out.setId(p1_in.getId());
                }


                void SingleCPUEnzymatic::perform(const particle_t &p1_in, const particle_t &p2_in, particle_t &p1_out, particle_t &p2_out) const {
                    if(p1_in.getType() == getCatalyst()) {
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

                void SingleCPUFission::perform(const particle_t &p1_in, const particle_t &p2_in, particle_t &p1_out, particle_t &p2_out) const {
                    // as long as the orientation is uniform, it does not matter of which type p1_in and p2_in are.
                    auto n3 = rand->getNormal3();
                    n3 /= sqrt(n3*n3);
                    p1_out.setType(getTo1());
                    p1_out.setPos(p1_in.getPos() + getWeight1()*getEductDistance()*n3);
                    p2_out.setType(getTo2());
                    p2_out.setPos(p1_in.getPos() - getWeight2()*getEductDistance()*n3);
                }

                void SingleCPUFusion::perform(const particle_t &p1_in, const particle_t &p2_in, particle_t &p1_out, particle_t &p2_out) const {
                    p1_out.setType(getTo());
                    if(getFrom1() == p1_in.getType()) {
                        p1_out.setPos(p1_in.getPos() - getWeight1()*(p2_in.getPos() - p1_in.getPos()));
                    } else {
                        p1_out.setPos(p2_in.getPos() - getWeight1()*(p2_in.getPos() - p1_in.getPos()));
                    }
                }
            }
        }
    }
}



