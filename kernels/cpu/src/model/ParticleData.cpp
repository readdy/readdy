/**
 * << detailed description >>
 *
 * @file ParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include "readdy/kernel/cpu/model/ParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

//
// Entry Impl
//

ParticleData::Entry::Entry(ParticleData::Entry &&rhs)
        : pos(std::move(rhs.pos)), force(std::move(rhs.force)), type(std::move(rhs.type)),
          deactivated(std::move(rhs.deactivated)), displacement(std::move(rhs.displacement)){ }

ParticleData::Entry &ParticleData::Entry::operator=(ParticleData::Entry &&rhs) {
    pos = std::move(rhs.pos);
    force = std::move(rhs.force);
    type = std::move(rhs.type);
    displacement = std::move(rhs.displacement);
    deactivated = std::move(rhs.deactivated);
    return *this;
}

const ParticleData::particle_type::pos_type &ParticleData::Entry::position() const {
    return pos;
}

bool ParticleData::Entry::is_deactivated() const {
    return deactivated;
}

//
// ParticleData Impl
//


ParticleData::ParticleData(readdy::model::KernelContext*const context)
        : blanks(std::vector<index_t>()), entries(), ids(), neighbors(), fixPos(context->getFixPositionFun())  { }

std::size_t ParticleData::size() const {
    return entries.size();
}

bool ParticleData::empty() const {
    return size() == getNDeactivated();
}

void ParticleData::clear() {
    while(!blanks.empty()) blanks.pop();
    entries.clear();
    ids.clear();
    neighbors.clear();
}

void ParticleData::addParticle(const ParticleData::particle_type &particle) {
    addParticles({particle});
}

void ParticleData::addParticles(const std::vector<ParticleData::particle_type> &particles) {
    for(const auto& p : particles) {
        if(!blanks.empty()) {
            const auto idx = blanks.top();
            blanks.pop();
            entries.at(idx) = {p};
            ids.at(idx) = p.getId();
        } else {
            entries.push_back({p});
            ids.push_back(p.getId());
            neighbors.push_back({});
        }
    }
}

void ParticleData::removeParticle(const ParticleData::particle_type &particle) {
    auto it_entries = begin();
    auto it_ids = ids.begin();
    std::size_t idx = 0;
    for(; it_entries != end(); ++it_entries, ++it_ids, ++idx) {
        if(!it_entries->is_deactivated() && *it_ids == particle.getId()) {
            blanks.push(idx);
            it_entries->deactivated = true;
            return;
        }
    }
    log::console()->error("Tried to remove particle ({}) which did not exist or was already deactivated!", particle);
}

void ParticleData::removeParticle(const ParticleData::index_t index) {
    auto& p = *(entries.begin() + index);
    if(!p.deactivated) {
        blanks.push(index);
        p.deactivated = true;
        // neighbors.at(index).clear();
    } else {
        log::console()->error("Tried to remove particle (index={}), that was already removed!", index);
    }
}

ParticleData::index_t ParticleData::getNDeactivated() const {
    return blanks.size();
}

readdy::model::Particle ParticleData::getParticle(const index_t index) const {
    const auto& entry = *(entries.begin() + index);
    if(entry.deactivated) {
        log::console()->error("Requested deactivated particle at index {}!", index);
    }
    return toParticle(entry, ids.at(index));
}

readdy::model::Particle ParticleData::toParticle(const Entry &e, const particle_type::id_type id) const {
    return readdy::model::Particle(e.pos, e.type, id);
}

ParticleData::index_t ParticleData::addEntry(ParticleData::EntryUpdate &&entry) {
    if(!blanks.empty()) {
        const auto idx = blanks.top();
        blanks.pop();
        ids.at(idx) = entry.id;
        entries.at(idx) = std::move(entry);
        neighbors.at(idx).clear();
        return idx;
    } else {
        ids.push_back(entry.id);
        entries.push_back(std::move(entry));
        neighbors.push_back({});
        return entries.size()-1;
    }
}

void ParticleData::removeEntry(index_t idx) {
    auto &entry = entries.at(idx);
    if(!entry.is_deactivated()) {
        entry.deactivated = true;
        blanks.push(idx);
    }
}

std::vector<ParticleData::index_t> ParticleData::update(update_t &&update_data) {
    std::vector<index_t> result;

    auto &&newEntries = std::move(std::get<0>(update_data));
    auto &&removedEntries = std::move(std::get<1>(update_data));
    result.reserve(newEntries.size());

    auto it_del = removedEntries.begin();
    for(auto&& newEntry : newEntries) {
        if(it_del != removedEntries.end()) {
            ids.at(*it_del) = newEntry.id;
            entries.at(*it_del) = std::move(newEntry);
            neighbors.at(*it_del).clear();
            result.push_back(*it_del);
            ++it_del;
        } else {
            result.push_back(addEntry(std::move(newEntry)));
        }
    }
    while(it_del != removedEntries.end()) {
        removeEntry(*it_del);
        ++it_del;
    }

    return result;
}

void ParticleData::setFixPosFun(const ctx_t::fix_pos_fun &f) {
    fixPos = f;
}

const ParticleData::particle_type::pos_type &ParticleData::pos(ParticleData::index_t idx) const {
    return entries.at(idx).pos;
}

void ParticleData::displace(ParticleData::Entry &entry, const readdy::model::Particle::pos_type &delta) {
    entry.pos += delta;
    fixPos(entry.pos);
    entry.displacement += sqrt(delta * delta);
}

ParticleData::Entry &ParticleData::entry_at(ParticleData::index_t idx) {
    return entries.at(idx);
}

const ParticleData::Entry &ParticleData::entry_at(ParticleData::index_t idx) const {
    return centry_at(idx);
}

const ParticleData::Entry &ParticleData::centry_at(ParticleData::index_t idx) const {
    return entries.at(idx);
}

ParticleData::neighbors_t &ParticleData::neighbors_at(ParticleData::index_t idx) {
    return neighbors.at(idx);
}

const ParticleData::neighbors_t &ParticleData::neighbors_at(ParticleData::index_t idx) const {
    return cneighbors_at(idx);
}

const ParticleData::neighbors_t &ParticleData::cneighbors_at(ParticleData::index_t idx) const {
    return neighbors.at(idx);
}


ParticleData::~ParticleData() = default;


ParticleData::EntryUpdate::EntryUpdate(const ParticleData::particle_type &particle) 
        : Entry(particle), id(particle.getId()) {}
}
}
}
}