/**
 * << detailed description >>
 *
 * @file ParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <readdy/kernel/cpu/model/ParticleData.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

ParticleData::ParticleData() : blanks(std::deque<index_t>()){
}

std::size_t ParticleData::size() const {
    // todo
    return 0; //ids.size() - blanks.size();
}

bool ParticleData::empty() const {
    return size() == 0;
}

void ParticleData::clear() {
    while(!blanks.empty()) blanks.pop();
}

void ParticleData::addParticle(const ParticleData::particle_type &particle) {
    addParticles({particle});
}

void ParticleData::addParticles(const std::vector<ParticleData::particle_type> &particles) {
    for(const auto& p : particles) {
        if(!blanks.empty()) {
            const auto idx = blanks.top();
            blanks.pop();
            //todo
        } else {
            //todo
        }
    }
}

ParticleData::~ParticleData() = default;


}
}
}
}