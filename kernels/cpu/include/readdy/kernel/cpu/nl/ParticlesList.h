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
 * @file ParticlesList.h
 * @brief << brief description >>
 * @author clonker
 * @date 24.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <vector>
#include <readdy/kernel/cpu/model/CPUParticleData.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class ParticlesList {
public:
    using particle_index = model::CPUParticleData::index_t;
    using particle_indices = std::vector<particle_index>;
    using mutex_type = std::mutex;
    using particles_lock = std::unique_lock<mutex_type>;

    ParticlesList() = default;

    ParticlesList(ParticlesList &&rhs) : _particles_mutex() {
        particles_lock rhs_lock(rhs._particles_mutex);
        _particles = std::move(rhs._particles);
    }

    ParticlesList &operator=(ParticlesList &&) = delete;

    ParticlesList(const ParticlesList &) = delete;

    ParticlesList &operator=(const ParticlesList &) = delete;

    virtual ~ParticlesList() = default;

    void add(const particle_index index) const {
        particles_lock lock(_particles_mutex);
        _particles.push_back(index);
    }

    void add(const particle_index index) {
        _particles.push_back(index);
    }

    void clear() {
        _particles.clear();
    }

    const particle_indices &get() const {
        return _particles;
    }

private:
    mutable std::mutex _particles_mutex{};
    mutable particle_indices _particles{};
};

}
}
}
}