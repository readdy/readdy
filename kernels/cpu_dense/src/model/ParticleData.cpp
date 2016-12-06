/**
 * << detailed description >>
 *
 * @file ParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */

#include <algorithm>
#include <readdy/common/logging.h>

#include "readdy/kernel/cpu_dense/model/SCPUParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace model {

ParticleData::ParticleData(const readdy::model::KernelContext *const ctx) : ParticleData(ctx, 0) {}

ParticleData::ParticleData(const readdy::model::KernelContext *const ctx, unsigned int capacity)
        : entries{}, n_marked(0), n_deactivated(0), deactivated_index(0), ctx(ctx) {
    entries.reserve(capacity);
}

ParticleData::iterator ParticleData::begin() {
    return entries.begin();
}

ParticleData::const_iterator ParticleData::begin() const {
    return cbegin();
}

ParticleData::const_iterator ParticleData::cbegin() const {
    return entries.cbegin();
}

ParticleData::iterator ParticleData::end() {
    return entries.begin() + deactivated_index;
}

ParticleData::const_iterator ParticleData::end() const {
    return cend();
}

ParticleData::const_iterator ParticleData::cend() const {
    return entries.cbegin() + deactivated_index;
}

std::size_t ParticleData::size() const {
    const auto s = n_marked.load();
    return s <= deactivated_index ? deactivated_index - s : 0;
}

std::size_t ParticleData::max_size() const {
    return entries.max_size();
}

bool ParticleData::empty() const {
    return size() == 0;
}

void ParticleData::clear() {
    deactivated_index = 0;
    n_deactivated = entries.size();
}

void ParticleData::addParticle(const ParticleData::particle_type &particle) {
    addParticles({particle});
}

void ParticleData::addParticles(const std::vector<ParticleData::particle_type> &particles) {
    auto added = particles.cbegin();
    auto it = entries.begin() + deactivated_index;
    while (added != particles.cend()) {
        if (n_deactivated > 0) {
            *it = {*added};

            --n_deactivated;
            ++deactivated_index;

            ++it;
        } else {
            entries.push_back({*added});
            ++deactivated_index;
        }
        ++added;
    }
}

void ParticleData::removeParticle(const ParticleData::particle_type &particle) {
    auto find_it = std::find_if(entries.begin(), entries.end(), [&](const Entry &e) {
        return !e.deactivated && e.id == particle.getId();
    });
    if (find_it != entries.end()) {
        removeParticle(static_cast<std::size_t>(std::distance(entries.begin(), find_it)));
    } else {
        log::console()->warn("Could not find and thus remove particle with id {}", particle.getId());
    }
}

void ParticleData::removeParticle(const std::size_t index) {
    entries.at(index).deactivated = true;

    std::swap(entries[index], entries[deactivated_index - 1]);

    ++n_deactivated;
    if (deactivated_index == 0) throw std::runtime_error("hier sollte man aber nicht hinkommen!1");
    --deactivated_index;
}

bool ParticleData::isMarkedForDeactivation(const std::size_t index) {
    return entries.at(index).deactivated;
}

std::size_t ParticleData::getDeactivatedIndex() const {
    return deactivated_index;
}

std::size_t ParticleData::getNDeactivated() const {
    return n_deactivated;
}

void ParticleData::deactivateMarked() {
    if (n_marked == 0) return;
    // sanity check: the deactivated_index is pointing to the
    // first (real) deactivated particle, i.e., marks the end of the
    // active data structure. "deactivated" is a vector<bool>
    // that is as long as the data, thus the deactivated_index
    // can be at most deactivated->begin() - deactivated->end().
    if (entries.size() < deactivated_index - 1) {
        throw std::runtime_error("this should not happen");
    }
    // we now are going backwards through the active part of the data structure,
    // starting with the first _active_ (but possible marked) particle
    auto deactivatedIt = begin() + deactivated_index - 1;
    // for each index in the markedForDeactivation data structure
    // (which is a set and thus sorted)
    for (auto it = begin(); it < end(); ++it) {
        if (it->deactivated) {
            const auto idx = it - begin();
            // if there are marked particles at the very end,
            // just shift the deactivated_index and increase n_deactivated
            while (deactivatedIt->deactivated && deactivatedIt != begin()) {
                --deactivated_index;
                ++n_deactivated;
                --deactivatedIt;
            }
            // since the deactivated_index might have decreased
            // so that we already have deactivated "idx", we check
            // if it has been deactivated already (by the above loop)
            if (idx < deactivated_index) {
                // performs swapping of this particle with the last active
                // particle
                removeParticle(static_cast<std::size_t>(idx));
                // if we are not at the begin already,
                // we want to decrease the current particle considered in
                // deactivatedIt
                if (deactivatedIt != begin()) --deactivatedIt;
            } else {
                // since the set is sorted and we start with the smallest idx,
                // we can stop here
                break;
            }
        }
    }
    n_marked = 0;
}

ParticleData::Entry &ParticleData::entry_at(const ParticleData::index_t idx) {
    return entries.at(idx);
}

const ParticleData::Entry &ParticleData::entry_at(const ParticleData::index_t idx) const {
    return entries.at(idx);
}

const ParticleData::Entry &ParticleData::centry_at(const ParticleData::index_t idx) const {
    return entries.at(idx);
}

readdy::model::Particle ParticleData::toParticle(const ParticleData::Entry &e) const {
    return {e.pos, e.type, e.id};
}

void ParticleData::update(ParticleData::update_t &&updates) {
    const auto& fixPos = ctx->getFixPositionFun();
    auto added = updates.begin();
    auto it = entries.begin() + deactivated_index;
    while (added != updates.end()) {
        fixPos((*added).pos);
        if (n_deactivated > 0) {
            *it = std::move(*added);

            --n_deactivated;
            ++deactivated_index;

            ++it;
        } else {
            entries.push_back(std::move(*added));
            ++deactivated_index;
        }
        ++added;
    }
}

void ParticleData::displace(ParticleData::Entry &entry, const readdy::model::Particle::pos_type &delta) {
    entry.pos += delta;
    ctx->getFixPositionFun()(entry.pos);
}

void ParticleData::deactivate(ParticleData::Entry &e) {
    e.deactivated = true;
    ++n_marked;
}

void ParticleData::deactivate(ParticleData::index_t idx) {
    deactivate(entry_at(idx));
}

ParticleData::~ParticleData() = default;

}
}
}
}