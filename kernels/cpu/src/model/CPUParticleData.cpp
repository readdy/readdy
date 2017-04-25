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
 * @file CPUParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <readdy/common/thread/Config.h>
#include <readdy/kernel/cpu/util/hilbert.h>
#include <readdy/kernel/cpu/util/config.h>
#include <readdy/common/Utils.h>
#include <readdy/common/signals.h>
#include "readdy/kernel/cpu/model/CPUParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

//
// Entry Impl
//

CPUParticleData::Entry::Entry(CPUParticleData::Entry &&rhs)
        : pos(std::move(rhs.pos)), force(std::move(rhs.force)), type(std::move(rhs.type)), id(std::move(rhs.id)),
          deactivated(std::move(rhs.deactivated)), displacement(std::move(rhs.displacement)){ }

CPUParticleData::Entry &CPUParticleData::Entry::operator=(CPUParticleData::Entry &&rhs) {
    pos = std::move(rhs.pos);
    force = std::move(rhs.force);
    type = std::move(rhs.type);
    id = std::move(rhs.id);
    displacement = std::move(rhs.displacement);
    deactivated = std::move(rhs.deactivated);
    return *this;
}

const CPUParticleData::particle_type::pos_type &CPUParticleData::Entry::position() const {
    return pos;
}

bool CPUParticleData::Entry::is_deactivated() const {
    return deactivated;
}

//
// CPUCPUParticleData Impl
//


CPUParticleData::CPUParticleData(readdy::model::KernelContext*const context, const readdy::util::thread::Config &_config)
        : blanks(std::vector<index_t>()), entries(), neighbors(), reorderSignal(std::make_unique<reorder_signal_t>()),
          fixPos(context->getFixPositionFun()), _config(_config), context(context) {}

std::size_t CPUParticleData::size() const {
    return entries.size();
}

bool CPUParticleData::empty() const {
    return size() == getNDeactivated();
}

void CPUParticleData::clear() {
    blanks.clear();
    entries.clear();
    neighbors.clear();
}

void CPUParticleData::addParticle(const CPUParticleData::particle_type &particle) {
    addParticles({particle});
}

void CPUParticleData::addParticles(const std::vector<CPUParticleData::particle_type> &particles) {
    for(const auto& p : particles) {
        if(!blanks.empty()) {
            const auto idx = blanks.back();
            blanks.pop_back();
            entries.at(idx) = {p};
        } else {
            entries.push_back({p});
            neighbors.push_back({});
        }
    }
}

void CPUParticleData::removeParticle(const CPUParticleData::particle_type &particle) {
    auto it_entries = begin();
    std::size_t idx = 0;
    for(; it_entries != end(); ++it_entries, ++idx) {
        if(!it_entries->is_deactivated() && it_entries->id == particle.getId()) {
            blanks.push_back(idx);
            it_entries->deactivated = true;
            return;
        }
    }
    log::error("Tried to remove particle ({}) which did not exist or was already deactivated!", particle);
}

void CPUParticleData::removeParticle(const CPUParticleData::index_t index) {
    auto& p = *(entries.begin() + index);
    if(!p.deactivated) {
        blanks.push_back(index);
        p.deactivated = true;
        // neighbors.at(index).clear();
    } else {
        log::error("Tried to remove particle (index={}), that was already removed!", index);
    }
}

CPUParticleData::index_t CPUParticleData::getNDeactivated() const {
    return blanks.size();
}

readdy::model::Particle CPUParticleData::getParticle(const index_t index) const {
    const auto& entry = *(entries.begin() + index);
    if(entry.deactivated) {
        log::error("Requested deactivated particle at index {}!", index);
    }
    return toParticle(entry);
}

readdy::model::Particle CPUParticleData::toParticle(const Entry &e) const {
    return readdy::model::Particle(e.pos, e.type, e.id);
}

CPUParticleData::index_t CPUParticleData::addEntry(CPUParticleData::Entry &&entry) {
    if(!blanks.empty()) {
        const auto idx = blanks.back();
        blanks.pop_back();
        entries.at(idx) = std::move(entry);
        neighbors.at(idx).clear();
        return idx;
    } else {
        entries.push_back(std::move(entry));
        neighbors.push_back({});
        return entries.size()-1;
    }
}

void CPUParticleData::removeEntry(index_t idx) {
    auto &entry = entries.at(idx);
    if(!entry.is_deactivated()) {
        entry.deactivated = true;
        blanks.push_back(idx);
    } else {
        log::critical("mist");
    }
}

std::vector<CPUParticleData::index_t> CPUParticleData::update(update_t &&update_data) {
    std::vector<index_t> result;

    auto &&newEntries = std::move(std::get<0>(update_data));
    auto &&removedEntries = std::move(std::get<1>(update_data));
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

void CPUParticleData::setFixPosFun(const ctx_t::fix_pos_fun &f) {
    fixPos = f;
}

const CPUParticleData::particle_type::pos_type &CPUParticleData::pos(CPUParticleData::index_t idx) const {
    return entries.at(idx).pos;
}

void CPUParticleData::displace(CPUParticleData::Entry &entry, const readdy::model::Particle::pos_type &delta) {
    entry.pos += delta;
    fixPos(entry.pos);
    entry.displacement += std::sqrt(delta * delta);
}

void CPUParticleData::blanks_moved_to_end() {
    auto n_blanks = blanks.size();
    std::iota(blanks.begin(), blanks.end(), size() - n_blanks - 1);
}

CPUParticleData::Entry &CPUParticleData::entry_at(CPUParticleData::index_t idx) {
    return entries.at(idx);
}

const CPUParticleData::Entry &CPUParticleData::entry_at(CPUParticleData::index_t idx) const {
    return centry_at(idx);
}

const CPUParticleData::Entry &CPUParticleData::centry_at(CPUParticleData::index_t idx) const {
    return entries.at(idx);
}

CPUParticleData::neighbors_t &CPUParticleData::neighbors_at(CPUParticleData::index_t idx) {
    return neighbors.at(idx);
}

const CPUParticleData::neighbors_t &CPUParticleData::neighbors_at(CPUParticleData::index_t idx) const {
    return cneighbors_at(idx);
}

const CPUParticleData::neighbors_t &CPUParticleData::cneighbors_at(CPUParticleData::index_t idx) const {
    return neighbors.at(idx);
}


CPUParticleData::~CPUParticleData() = default;


CPUParticleData::iterator CPUParticleData::begin() {
    return entries.begin();
}

CPUParticleData::iterator CPUParticleData::end() {
    return entries.end();
}

CPUParticleData::const_iterator CPUParticleData::cbegin() const {
    return entries.cbegin();
}

CPUParticleData::const_iterator CPUParticleData::cend() const {
    return entries.cend();
}


CPUParticleData::const_iterator CPUParticleData::begin() const {
    return entries.cbegin();
}


CPUParticleData::const_iterator CPUParticleData::end() const {
    return entries.cend();
}

void CPUParticleData::blanks_moved_to_front() {
    std::iota(blanks.begin(), blanks.end(), 0);
}

CPUParticleData::index_t CPUParticleData::getIndexForId(const particle_type::id_type id) const {
    auto find_it = std::find_if(entries.begin(), entries.end(), [id](const Entry& e) {
        return !e.is_deactivated() && e.id == id;
    });
    if(find_it != entries.end()) {
        return static_cast<index_t>(std::distance(entries.begin(), find_it));
    }
    throw std::out_of_range("requested id was not to be found in particle data");
}

std::vector<std::size_t>
CPUParticleData::addTopologyParticles(const std::vector<CPUParticleData::top_particle_type> &particles) {
    std::vector<entries_t::size_type> indices;
    indices.reserve(particles.size());
    for(const auto& p : particles) {
        if(!blanks.empty()) {
            const auto idx = blanks.back();
            blanks.pop_back();
            entries.at(idx) = {p};
            indices.push_back(idx);
        } else {
            entries.push_back({p});
            neighbors.push_back({});
            indices.push_back(entries.size()-1);
        }
    }
    return indices;
}

void CPUParticleData::reserve(std::size_t n) {
    entries.reserve(n);
    neighbors.reserve(n);
}

const CPUParticleData::ctx_t::fix_pos_fun &CPUParticleData::fixPosFun() const {
    return fixPos;
}

void CPUParticleData::setPBCFun(const CPUParticleData::ctx_t::pbc_fun &f) {
    pbc = f;
}

const CPUParticleData::ctx_t::pbc_fun &CPUParticleData::pbcFun() const {
    return pbc;
}

void CPUParticleData::hilbert_sort(const scalar grid_width) {
    if(!empty()) {
        using indices_it = std::vector<std::size_t>::iterator;
        std::vector<std::size_t> hilbert_indices;
        std::vector<std::size_t> indices(size());

        const auto grainSize = size() / _config.nThreads();
        std::iota(indices.begin(), indices.end(), 0);
        hilbert_indices.resize(size());

        auto worker = [&](const_iterator begin, const_iterator end, indices_it hilbert_begin) {
            {
                auto it = begin;
                auto hilbert_it = hilbert_begin;
                for (; it != end; ++it, ++hilbert_it) {
                    if (!it->is_deactivated()) {
                        const auto integer_coordinates = project(it->position(), grid_width);
                        *hilbert_it = 1 + static_cast<std::size_t>(hilbert_c2i(3, CHAR_BIT, integer_coordinates.data()));
                    } else {
                        *hilbert_it = 0;
                    }
                }
            }
        };
        {
            std::vector<threading_model> threads;
            threads.reserve(_config.nThreads());
            auto data_it = begin();
            auto hilberts_it = hilbert_indices.begin();
            for (std::size_t i = 0; i < _config.nThreads() - 1; ++i) {
                threads.emplace_back(worker, data_it, data_it + grainSize, hilberts_it);
                data_it += grainSize;
                hilberts_it += grainSize;
            }
            threads.emplace_back(worker, data_it, end(), hilberts_it);
        }
        {
            std::sort(indices.begin(), indices.end(),
                      [&hilbert_indices](std::size_t i, std::size_t j) -> bool {
                          return hilbert_indices[i] < hilbert_indices[j];
                      });
        }

        blanks_moved_to_front();
        std::vector<std::size_t> inverseIndices(indices.size());
        for(std::size_t i = 0; i < indices.size(); ++i) {
            inverseIndices[indices[i]] = i;
        }
        reorderSignal->fire_signal(inverseIndices);
        readdy::util::collections::reorder_destructive(inverseIndices.begin(), inverseIndices.end(), begin());
    }
}

std::array<unsigned long long, 3> CPUParticleData::project(vec3 value, scalar grid_width) const{
    const auto& box_size = context->getBoxSize();
    const auto i = static_cast<const unsigned long long>((value.x + .5 * box_size[0]) / grid_width);
    const auto j = static_cast<const unsigned long long>((value.y + .5 * box_size[1]) / grid_width);
    const auto k = static_cast<const unsigned long long>((value.z + .5 * box_size[2]) / grid_width);
    return {{i, j, k}};
}

readdy::signals::scoped_connection
CPUParticleData::registerReorderEventListener(const reorder_signal_t::slot_type &slot) {
    return reorderSignal->connect_scoped(slot);
}

}
}
}
}