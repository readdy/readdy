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
 * @file SingleCPUParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#include <numeric>
#include <readdy/common/make_unique.h>
#include <readdy/common/logging.h>
#include <readdy/kernel/singlecpu/model/SCPUParticleData.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {

SCPUParticleData::SCPUParticleData(bool useMarkedSet)
        : SCPUParticleData(0, useMarkedSet) {}

SCPUParticleData::SCPUParticleData(unsigned int capacity)
        : SCPUParticleData(capacity, true) {}

SCPUParticleData::SCPUParticleData() : SCPUParticleData(0, true) {
}

SCPUParticleData::SCPUParticleData(unsigned int capacity, bool useMarkedSet)
        : useMarkedSet(useMarkedSet) {
    ids = std::make_unique<std::vector<readdy::model::Particle::id_type>>(capacity);
    positions = std::make_unique<std::vector<readdy::model::Vec3>>(capacity);
    forces = std::make_unique<std::vector<readdy::model::Vec3>>(capacity);
    type = std::make_unique<std::vector<unsigned int>>(capacity);
    if (useMarkedSet) {
        markedForDeactivation = std::make_unique<std::set<size_t>>();
    } else {
        n_marked = 0;
    }
    deactivated = std::make_unique<std::vector<char>>(capacity);
    std::fill(deactivated->begin(), deactivated->end(), true);
    n_deactivated = capacity;
    deactivated_index = 0;
}

void SCPUParticleData::swap(SCPUParticleData &rhs) {
    std::swap(ids, rhs.ids);
    std::swap(positions, rhs.positions);
    std::swap(forces, rhs.forces);
    std::swap(type, rhs.type);
    std::swap(deactivated, rhs.deactivated);
    std::swap(deactivated_index, rhs.deactivated_index);
    std::swap(n_deactivated, rhs.n_deactivated);
    if (useMarkedSet) {
        std::swap(markedForDeactivation, rhs.markedForDeactivation);
    } else {
        n_marked = rhs.n_marked.exchange(n_marked);
    }
    std::swap(useMarkedSet, rhs.useMarkedSet);
}

size_t SCPUParticleData::size() const {
    auto s = useMarkedSet ? markedForDeactivation->size() : n_marked.load();
    return s <= deactivated_index ? deactivated_index - s : 0;
}

size_t SCPUParticleData::max_size() const {
    return ids->max_size();
}

bool SCPUParticleData::empty() const {
    return size() == 0;
}

void SCPUParticleData::addParticle(const readdy::model::Particle &particle) {
    addParticles({particle});
};

void SCPUParticleData::addParticles(const std::vector<readdy::model::Particle> &particles) {
    auto added = particles.cbegin();
    auto ids_it = ids->begin() + deactivated_index;
    auto positions_it = positions->begin() + deactivated_index;
    auto forces_it = forces->begin() + deactivated_index;
    auto type_it = type->begin() + deactivated_index;
    auto deactivated_it = deactivated->begin() + deactivated_index;
    while (added != particles.cend()) {
        if (n_deactivated > 0) {

            *ids_it = added->getId();
            *positions_it = added->getPos();
            *forces_it = {0, 0, 0};
            *type_it = added->getType();
            *deactivated_it = false;

            --n_deactivated;
            ++deactivated_index;

            ++ids_it;
            ++positions_it;
            ++forces_it;
            ++type_it;
            ++deactivated_it;
        } else {
            ids->push_back(added->getId());
            positions->push_back(added->getPos());
            forces->push_back({0, 0, 0});
            type->push_back(added->getType());
            deactivated->push_back(false);
            ++deactivated_index;
        }
        ++added;
    }
}

void SCPUParticleData::removeParticle(const size_t index) {
    (*deactivated)[index] = true;

    std::swap((*ids)[index], (*ids)[deactivated_index - 1]);
    std::swap((*positions)[index], (*positions)[deactivated_index - 1]);
    std::swap((*forces)[index], (*forces)[deactivated_index - 1]);
    std::swap((*type)[index], (*type)[deactivated_index - 1]);
    std::swap((*deactivated)[index], (*deactivated)[deactivated_index - 1]);

    ++n_deactivated;
    if (deactivated_index == 0) throw std::runtime_error("hier sollte man aber nicht hinkommen!1");
    --deactivated_index;
}

void SCPUParticleData::markForDeactivation(size_t index) {
    auto it = begin_deactivated() + index;
    if (useMarkedSet) {
        std::lock_guard<std::mutex> lock(markedForDeactivationMutex);
        markedForDeactivation->insert(index);
    } else {
        if (*it == false) {
            ++n_marked;
        } else {
            log::console()->error("this should not have happened! (idx={})", index);
        }
    }
    *it = true;
}

void SCPUParticleData::deactivateMarked() {
    if (useMarkedSet) {
        deactivateMarkedSet();
    } else {
        deactivateMarkedNoSet();
    }
}

void SCPUParticleData::removeParticle(const readdy::model::Particle &particle) {
    auto &&beginIt = begin_ids();
    auto &&endIt = end_ids();
    auto &&it = std::find(beginIt, endIt, particle.getId());
    if (it != endIt) {
        removeParticle(it - beginIt);
    } else {
        log::console()->warn("Could not find and thus remove particle");
    }
}

std::vector<readdy::model::Particle::id_type>::iterator SCPUParticleData::begin_ids() {
    return ids->begin();
}

std::vector<readdy::model::Particle::id_type>::iterator SCPUParticleData::end_ids() {
    return ids->begin() + deactivated_index;
}

std::vector<readdy::model::Vec3>::iterator SCPUParticleData::begin_positions() {
    return positions->begin();
}

std::vector<readdy::model::Vec3>::iterator SCPUParticleData::end_positions() {
    return positions->begin() + deactivated_index;
}

std::vector<readdy::model::Vec3>::iterator SCPUParticleData::begin_forces() {
    return forces->begin();
}

std::vector<readdy::model::Vec3>::iterator SCPUParticleData::end_forces() {
    return forces->begin() + deactivated_index;
}

std::vector<unsigned int>::iterator SCPUParticleData::begin_types() {
    return type->begin();
}

std::vector<unsigned int>::iterator SCPUParticleData::end_types() {
    return type->begin() + deactivated_index;
}

readdy::model::Particle SCPUParticleData::operator[](const size_t index) const {
    return readdy::model::Particle(*(begin_positions() + index), *(begin_types() + index), *(begin_ids() + index));
}


SCPUParticleData &SCPUParticleData::operator=(SCPUParticleData &&rhs) {
    if (this != &rhs) {
        std::unique_lock<std::mutex> lhs_lock(markedForDeactivationMutex, std::defer_lock);
        std::unique_lock<std::mutex> rhs_lock(rhs.markedForDeactivationMutex, std::defer_lock);
        std::lock(lhs_lock, rhs_lock);
        ids = std::move(rhs.ids);
        positions = std::move(rhs.positions);
        forces = std::move(rhs.forces);
        type = std::move(rhs.type);
        deactivated = std::move(rhs.deactivated);
        deactivated_index = std::move(rhs.deactivated_index);
        n_deactivated = std::move(rhs.n_deactivated);
        if (useMarkedSet) {
            markedForDeactivation = std::move(rhs.markedForDeactivation);
        } else {
            n_marked = rhs.n_marked.load();
        }
        useMarkedSet = std::move(rhs.useMarkedSet);
    }
    return *this;
};

SCPUParticleData::SCPUParticleData(SCPUParticleData &&rhs) {
    std::unique_lock<std::mutex> rhs_lock(rhs.markedForDeactivationMutex);
    ids = std::move(rhs.ids);
    positions = std::move(rhs.positions);
    forces = std::move(rhs.forces);
    type = std::move(rhs.type);
    deactivated = std::move(rhs.deactivated);
    deactivated_index = std::move(rhs.deactivated_index);
    n_deactivated = std::move(rhs.n_deactivated);
    if (useMarkedSet) {
        markedForDeactivation = std::move(rhs.markedForDeactivation);
    } else {
        n_marked = rhs.n_marked.load();
    }
    useMarkedSet = std::move(rhs.useMarkedSet);
};

SCPUParticleData::~SCPUParticleData() {
}

std::vector<readdy::model::Particle::id_type>::const_iterator SCPUParticleData::begin_ids() const {
    return cbegin_ids();
}

std::vector<readdy::model::Particle::id_type>::const_iterator SCPUParticleData::cbegin_ids() const {
    return ids->cbegin();
}

std::vector<readdy::model::Particle::id_type>::const_iterator SCPUParticleData::end_ids() const {
    return cend_ids();
}

std::vector<readdy::model::Particle::id_type>::const_iterator SCPUParticleData::cend_ids() const {
    return ids->begin() + deactivated_index;
}

std::vector<readdy::model::Vec3>::const_iterator SCPUParticleData::begin_positions() const {
    return cbegin_positions();
}

std::vector<readdy::model::Vec3>::const_iterator SCPUParticleData::cbegin_positions() const {
    return positions->cbegin();
}

std::vector<readdy::model::Vec3>::const_iterator SCPUParticleData::end_positions() const {
    return cend_positions();
}

std::vector<readdy::model::Vec3>::const_iterator SCPUParticleData::cend_positions() const {
    return positions->cbegin() + deactivated_index;
}

std::vector<readdy::model::Vec3>::const_iterator SCPUParticleData::begin_forces() const {
    return cbegin_forces();
}

std::vector<readdy::model::Vec3>::const_iterator SCPUParticleData::cbegin_forces() const {
    return forces->cbegin();
}

std::vector<readdy::model::Vec3>::const_iterator SCPUParticleData::end_forces() const {
    return cend_forces();
}

std::vector<readdy::model::Vec3>::const_iterator SCPUParticleData::cend_forces() const {
    return forces->cbegin() + deactivated_index;
}

std::vector<unsigned int>::const_iterator SCPUParticleData::begin_types() const {
    return cbegin_types();
}

std::vector<unsigned int>::const_iterator SCPUParticleData::cbegin_types() const {
    return type->cbegin();
}

std::vector<unsigned int>::const_iterator SCPUParticleData::end_types() const {
    return cend_types();
}

std::vector<unsigned int>::const_iterator SCPUParticleData::cend_types() const {
    return type->cbegin() + deactivated_index;
}

bool SCPUParticleData::isMarkedForDeactivation(const size_t index) {
    return (*deactivated)[index];
}

size_t SCPUParticleData::getDeactivatedIndex() const {
    return deactivated_index;
}

size_t SCPUParticleData::getNDeactivated() const {
    return n_deactivated;
}

void SCPUParticleData::setParticleData(const readdy::model::Particle &particle, const size_t &index) {
    (*ids)[index] = particle.getId();
    (*positions)[index] = particle.getPos();
    (*type)[index] = particle.getType();
}

void SCPUParticleData::clear() {
    deactivated_index = 0;
    n_deactivated = ids->size();
}

std::vector<char>::iterator SCPUParticleData::begin_deactivated() {
    return deactivated->begin();
}

std::vector<char>::const_iterator SCPUParticleData::begin_deactivated() const {
    return cbegin_deactivated();
}

std::vector<char>::iterator SCPUParticleData::end_deactivated() {
    return deactivated->begin() + deactivated_index;
}

std::vector<char>::const_iterator SCPUParticleData::end_deactivated() const {
    return cend_deactivated();
}

std::vector<char>::const_iterator SCPUParticleData::cend_deactivated() const {
    return deactivated->cbegin() + deactivated_index;
}

std::vector<char>::const_iterator SCPUParticleData::cbegin_deactivated() const {
    return deactivated->cbegin();
}

void SCPUParticleData::deactivateMarkedNoSet() {
    if (n_marked == 0) return;
    // sanity check: the deactivated_index is pointing to the
    // first (real) deactivated particle, i.e., marks the end of the
    // active data structure. "deactivated" is a vector<bool>
    // that is as long as the data, thus the deactivated_index
    // can be at most deactivated->begin() - deactivated->end().
    if (deactivated->size() < deactivated_index - 1) {
        throw std::runtime_error("this should not happen");
    }
    // we now are going backwards through the active part of the data structure,
    // starting with the first _active_ (but possible marked) particle
    auto deactivatedIt = deactivated->begin() + deactivated_index - 1;
    // for each index in the markedForDeactivation data structure
    // (which is a set and thus sorted)
    for (auto it = begin_deactivated(); it < end_deactivated(); ++it) {
        if (*it) {
            const auto idx = it - begin_deactivated();
            // if there are marked particles at the very end,
            // just shift the deactivated_index and increase n_deactivated
            while (*deactivatedIt && deactivatedIt != deactivated->begin()) {
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
                removeParticle(idx);
                // if we are not at the begin already,
                // we want to decrease the current particle considered in
                // deactivatedIt
                if (deactivatedIt != deactivated->begin()) --deactivatedIt;
            } else {
                // since the set is sorted and we start with the smallest idx,
                // we can stop here
                break;
            }
        }
    }
    n_marked = 0;
}

void SCPUParticleData::deactivateMarkedSet() {
    // if we havent marked anything, return
    if (markedForDeactivation->size() == 0) return;
    // sanity check: the deactivated_index is pointing to the
    // first (real) deactivated particle, i.e., marks the end of the
    // active data structure. "deactivated" is a vector<bool>
    // that is as long as the data, thus the deactivated_index
    // can be at most deactivated->begin() - deactivated->end().
    if (deactivated->size() < deactivated_index - 1) {
        throw std::runtime_error("this should not happen");
    }
    // if we have active particles
    if (deactivated_index > 0) {
        // we now are going backwards through the active part of the data structure,
        // starting with the first _active_ (but possible marked) particle
        auto deactivatedIt = deactivated->begin() + deactivated_index - 1;
        // for each index in the markedForDeactivation data structure
        // (which is a set and thus sorted)
        for (auto &&idx : *markedForDeactivation) {
            // if there are marked particles at the very end,
            // just shift the deactivated_index and increase n_deactivated
            while (*deactivatedIt && deactivatedIt != deactivated->begin()) {
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
                removeParticle(idx);
                // if we are not at the begin already,
                // we want to decrease the current particle considered in
                // deactivatedIt
                if (deactivatedIt != deactivated->begin()) --deactivatedIt;
            } else {
                // since the set is sorted and we start with the smallest idx,
                // we can stop here
                break;
            }
        }
    }
    markedForDeactivation->clear();
}
}
}
}
}




