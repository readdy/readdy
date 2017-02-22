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
 * @file ParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */

#include <algorithm>
#include <readdy/common/logging.h>

#include "readdy/kernel/cpu_dense/model/CPUDParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace model {

CPUDParticleData::CPUDParticleData(const readdy::model::KernelContext *const ctx) : CPUDParticleData(ctx, 0) {}

CPUDParticleData::CPUDParticleData(const readdy::model::KernelContext *const ctx, unsigned int capacity)
        : entries{}, n_marked(0), n_deactivated(0), deactivated_index(0), ctx(ctx) {
    entries.reserve(capacity);
}

CPUDParticleData::iterator CPUDParticleData::begin() {
    return entries.begin();
}

CPUDParticleData::const_iterator CPUDParticleData::begin() const {
    return cbegin();
}

CPUDParticleData::const_iterator CPUDParticleData::cbegin() const {
    return entries.cbegin();
}

CPUDParticleData::iterator CPUDParticleData::end() {
    return entries.begin() + deactivated_index;
}

CPUDParticleData::const_iterator CPUDParticleData::end() const {
    return cend();
}

CPUDParticleData::const_iterator CPUDParticleData::cend() const {
    return entries.cbegin() + deactivated_index;
}

std::size_t CPUDParticleData::size() const {
    const auto s = n_marked.load();
    return s <= deactivated_index ? deactivated_index - s : 0;
}

std::size_t CPUDParticleData::max_size() const {
    return entries.max_size();
}

bool CPUDParticleData::empty() const {
    return size() == 0;
}

void CPUDParticleData::clear() {
    deactivated_index = 0;
    n_deactivated = entries.size();
}

void CPUDParticleData::addParticle(const CPUDParticleData::particle_type &particle) {
    addParticles({particle});
}

void CPUDParticleData::addParticles(const std::vector<CPUDParticleData::particle_type> &particles) {
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

void CPUDParticleData::removeParticle(const CPUDParticleData::particle_type &particle) {
    auto find_it = std::find_if(entries.begin(), entries.end(), [&](const Entry &e) {
        return !e.deactivated && e.id == particle.getId();
    });
    if (find_it != entries.end()) {
        removeParticle(static_cast<std::size_t>(std::distance(entries.begin(), find_it)));
    } else {
        log::warn("Could not find and thus remove particle with id {}", particle.getId());
    }
}

void CPUDParticleData::removeParticle(const std::size_t index) {
    entries.at(index).deactivated = true;

    std::swap(entries[index], entries[deactivated_index - 1]);

    ++n_deactivated;
    if (deactivated_index == 0) throw std::runtime_error("hier sollte man aber nicht hinkommen!1");
    --deactivated_index;
}

bool CPUDParticleData::isMarkedForDeactivation(const std::size_t index) {
    return entries.at(index).deactivated;
}

std::size_t CPUDParticleData::getDeactivatedIndex() const {
    return deactivated_index;
}

std::size_t CPUDParticleData::getNDeactivated() const {
    return n_deactivated;
}

void CPUDParticleData::deactivateMarked() {
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

CPUDParticleData::Entry &CPUDParticleData::entry_at(const CPUDParticleData::index_t idx) {
    return entries.at(idx);
}

const CPUDParticleData::Entry &CPUDParticleData::entry_at(const CPUDParticleData::index_t idx) const {
    return entries.at(idx);
}

const CPUDParticleData::Entry &CPUDParticleData::centry_at(const CPUDParticleData::index_t idx) const {
    return entries.at(idx);
}

readdy::model::Particle CPUDParticleData::toParticle(const CPUDParticleData::Entry &e) const {
    return {e.pos, e.type, e.id};
}

void CPUDParticleData::update(CPUDParticleData::update_t &&updates) {
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

void CPUDParticleData::displace(CPUDParticleData::Entry &entry, const readdy::model::Particle::pos_type &delta) {
    entry.pos += delta;
    ctx->getFixPositionFun()(entry.pos);
}

void CPUDParticleData::deactivate(CPUDParticleData::Entry &e) {
    e.deactivated = true;
    ++n_marked;
}

void CPUDParticleData::deactivate(CPUDParticleData::index_t idx) {
    deactivate(entry_at(idx));
}

CPUDParticleData::~CPUDParticleData() = default;

}
}
}
}