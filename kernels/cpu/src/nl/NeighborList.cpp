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
 * @file NeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 24.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu/nl/NeighborList.h>
#include <readdy/kernel/cpu/util/config.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

NeighborList::NeighborList(const model::CPUParticleData &data, const readdy::model::KernelContext &context,
                           const readdy::util::thread::Config &config, NeighborList::skin_size_t skin)
        : _data(data), _context(context), _config(config), _skin(skin), _cell_container(data, context, config) {
}

const CellContainer &NeighborList::cell_container() const {
    return _cell_container;
}

void NeighborList::set_up() {
    _max_cutoff = calculate_max_cutoff();
    if(_max_cutoff > 0) {
        _cell_container.subdivide(_max_cutoff + _skin);
        _cell_container.refine_uniformly();
        _cell_container.setup_uniform_neighbors();
        fill_container();
    }
}

scalar NeighborList::calculate_max_cutoff() {
    scalar max_cutoff = 0;
    for(const auto& entry : _context.potentials().potentials_order2()) {
        for(const auto& potential : entry.second) {
            max_cutoff = max_cutoff < potential->getCutoffRadius() ? potential->getCutoffRadius() : max_cutoff;
        }
    }
    for(const auto& entry : _context.reactions().order2()) {
        for(const auto& reaction : entry.second) {
            max_cutoff = max_cutoff < reaction->getEductDistance() ? reaction->getEductDistance() : max_cutoff;
        }
    }
    return max_cutoff;
}

void NeighborList::update() {
    if(_max_cutoff > 0) {
        // todo update as needed
    }
}

const model::CPUParticleData &NeighborList::data() const {
    return _data;
}

const readdy::util::thread::Config &NeighborList::config() const {
    return _config;
}

void NeighborList::fill_container() {
    const auto grainSize = _data.size() / _config.nThreads();
    const auto data_begin = _data.cbegin();
    auto worker = [=](model::CPUParticleData::const_iterator begin, model::CPUParticleData::const_iterator end, const CellContainer& container) {
        for (auto it = begin; it != end; ++it) {
            if(!it->is_deactivated()) {
                container.insert_particle(static_cast<CellContainer::particle_index>(std::distance(data_begin, it)));
            }
        }
    };
    std::vector<threading_model> threads;
    threads.reserve(_config.nThreads());
    auto it = _data.cbegin();
    for (auto i = 0; i < _config.nThreads() - 1; ++i) {
        threads.emplace_back(worker, it, it + grainSize, std::cref(_cell_container));
        it += grainSize;
    }
    threads.emplace_back(worker, it, _data.end(), std::cref(_cell_container));
    /*std::size_t i = 0;
    for(auto it = _data.begin(); it != _data.end(); ++it, ++i) {
        if(!it->is_deactivated()) {
            _cell_container.insert_particle(i);
        }
    }*/
}


}
}
}
}
