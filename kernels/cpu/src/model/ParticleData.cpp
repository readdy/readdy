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

//
// Neighbor Impl
//

ParticleData::Neighbor::Neighbor(const index_t idx, const double d2) : idx(idx), d2(d2) {}

ParticleData::Neighbor::Neighbor(Neighbor &&rhs) : idx(rhs.idx), d2(std::move(rhs.d2)) {
}

ParticleData::Neighbor& ParticleData::Neighbor::operator=(ParticleData::Neighbor &&rhs) {
    idx = rhs.idx;
    d2 = std::move(rhs.d2);
    return *this;
}

//
// Entry Impl
//

ParticleData::Entry::Entry(ParticleData::Entry &&rhs)
        : id(std::move(rhs.id)), pos(std::move(rhs.pos)), force(std::move(rhs.force)), type(std::move(rhs.type)),
          deactivated(std::move(rhs.deactivated)), displacement(std::move(rhs.displacement)){ }

ParticleData::Entry &ParticleData::Entry::operator=(ParticleData::Entry &&rhs) {
    id = std::move(rhs.id);
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
        : blanks(std::vector<index_t>()), entries(), neighbors(), fixPos(context->getFixPositionFun())  { }

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
            entries[idx] = {p};
        } else {
            entries.push_back({p});
            neighbors.push_back({});
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

ParticleData::index_t ParticleData::addEntry(ParticleData::Entry &&entry) {
    if(!blanks.empty()) {
        const auto idx = blanks.top();
        blanks.pop();
        entries.at(idx) = std::move(entry);
        return idx;
    } else {
        entries.push_back(std::move(entry));
        neighbors.push_back({});
        return entries.size();
    }
}

void ParticleData::removeEntry(index_t idx) {
    auto& entry = entries.at(idx);
    if(!entry.is_deactivated()) {
        entry.deactivated = true;
        blanks.push(idx);
    }
}

ParticleData::index_t ParticleData::getEntryIndex(const ParticleData::Entry *const entry) const {
    return static_cast<index_t>(entry - &*entries.begin());
}

std::vector<ParticleData::index_t> ParticleData::update(update_t &&update_data) {
    auto &&newEntries = std::move(std::get<0>(update_data));
    auto &&removedEntries = std::move(std::get<1>(update_data));

    std::vector<index_t> result;
    result.reserve(newEntries.size());

    auto it_del = removedEntries.begin();
    for(auto&& newEntry : newEntries) {
        if(it_del != removedEntries.end()) {
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

void ParticleData::displace(index_t index, const particle_type::pos_type &delta) {
    displace(entries.at(index), delta);
}

void ParticleData::setPosition(index_t idx, particle_type::pos_type &&newPosition) {
    auto& entry = entries.at(idx);
    const auto delta = newPosition - entry.pos;
    entry.pos = std::move(newPosition);
    fixPos(entry.pos);
    entry.displacement += sqrt(delta * delta);
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


ParticleData::~ParticleData() = default;


}
}
}
}