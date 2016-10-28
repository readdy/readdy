/**
 * << detailed description >>
 *
 * @file ParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <readdy/kernel/cpu/model/ParticleData.h>
#include <readdy/common/logging.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

ParticleData::ParticleData() : blanks(std::deque<index_t>()), entries({}){
}

std::size_t ParticleData::size() const {
    return entries.size() - blanks.size();
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
            entries[idx] = p;
        } else {
            entries.push_back({p});
        }
    }
}

void ParticleData::removeParticle(const ParticleData::particle_type &particle) {
    auto it = std::find_if(entries.begin(), entries.end(), [&particle](const Entry& e) {
        return e.id == particle.getId() && !e.deactivated;
    });
    if(it != entries.end()) {
        blanks.push(it - entries.begin());
        (*it).deactivated = true;
    } else {
        log::console()->error("Tried to remove particle ({}) which did not exist or was already deactivated!", particle);
    }
}

void ParticleData::removeParticle(const ParticleData::index_t index) {
    auto& p = *(entries.begin() + index);
    if(!p.deactivated) {
        blanks.push(index);
        p.deactivated = true;
    } else {
        log::console()->error("Tried to remove particle (index={}), that was already removed!", index);
    }
}

ParticleData::index_t ParticleData::getNDeactivated() const {
    return blanks.size();
}

readdy::model::Particle ParticleData::getParticle(const ParticleData::index_t index) const {
    const auto& entry = *(entries.begin() + index);
    if(entry.deactivated) {
        log::console()->error("Requested deactivated particle!");
    }
    return toParticle(entry);
}

readdy::model::Particle ParticleData::toParticle(const ParticleData::Entry &e) const {
    return readdy::model::Particle(e.pos, e.type, e.id);
}

void ParticleData::addEntries(const std::vector<ParticleData::Entry> &newEntries) {
    for(const auto& entry : newEntries) {
        addEntry(entry);
    }
}

ParticleData::index_t ParticleData::addEntry(ParticleData::Entry entry) {
    if(!blanks.empty()) {
        const auto idx = blanks.top();
        blanks.pop();
        entries[idx] = std::move(entry);
        return idx;
    } else {
        entries.push_back(std::move(entry));
        return entries.size()-1;
    }
}

ParticleData::~ParticleData() = default;


}
}
}
}