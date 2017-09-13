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
 *
 *
 * @file CLLNeighborList.cpp
 * @brief 
 * @author clonker
 * @date 9/13/17
 */

#include <readdy/kernel/cpu/nl/CLLNeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {


ContiguousCLLNeighborList::ContiguousCLLNeighborList(model::CPUParticleData &data,
                                                     const readdy::model::KernelContext &context,
                                                     const readdy::util::thread::Config &config)
        : NeighborList(data, context, config), ccll(data, context, config) {}

void ContiguousCLLNeighborList::set_up(const util::PerformanceNode &node) {
    auto t = node.timeit();
    _max_cutoff = _context.get().calculateMaxCutoff();
    _max_cutoff_skin_squared = (_max_cutoff + _skin) * (_max_cutoff + _skin);

    ccll.setUp(_skin, cll_radius, node.subnode("setUp CLL"));
    fill_verlet_list(node.subnode("fill verlet list"));
    _is_set_up = true;
}

void ContiguousCLLNeighborList::fill_verlet_list(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if(_max_cutoff > 0) {
        auto cix = ccll.cellIndex();
        auto bix = ccll.binsIndex();
        auto nix = ccll.neighborIndex();

        const auto grainSize = cix.size() / _config.get().nThreads();
        auto worker = [this, cix, bix, nix](std::size_t tid, std::size_t begin, std::size_t end) {
            const auto &d2 = _context.get().distSquaredFun();
            auto &data = _data.get();

            for(std::size_t cellIndex = begin; cellIndex < end; ++cellIndex) {

                for(auto itParticles = ccll.particlesBegin(cellIndex); itParticles != ccll.particlesEnd(cellIndex); ++itParticles) {
                    auto particle = *itParticles;
                    auto &entry = data.entry_at(particle);
                    auto &neighbors = data.neighbors_at(particle);
                    neighbors.clear();

                    for(auto itPP = ccll.particlesBegin(cellIndex); itPP != ccll.particlesEnd(cellIndex); ++itPP) {
                        auto pparticle = *itPP;
                        if(particle != pparticle) {
                            const auto &pp = data.entry_at(pparticle);
                            if(!pp.deactivated) {
                                const auto distSquared = d2(entry.pos, pp.pos);
                                if (distSquared < _max_cutoff_skin_squared) {
                                    neighbors.push_back(pparticle);
                                }
                            }
                        }
                    }

                    for(auto itNeighborCell = ccll.neighborsBegin(cellIndex); itNeighborCell != ccll.neighborsEnd(cellIndex); ++itNeighborCell) {
                        for(auto itNeighborParticle = ccll.particlesBegin(*itNeighborCell); itNeighborParticle != ccll.particlesEnd(*itNeighborCell); ++itNeighborParticle) {
                            const auto &neighbor = data.entry_at(*itNeighborParticle);
                            if(!neighbor.deactivated) {
                                const auto distSquared = d2(entry.pos, neighbor.pos);
                                if (distSquared < _max_cutoff_skin_squared) {
                                    neighbors.push_back(*itNeighborParticle);
                                }
                            }
                        }
                    }
                }
            }
        };
        const auto& executor = *_config.get().executor();
        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(_config.get().nThreads());
        auto it = 0_z;
        for(int i = 0; i < _config.get().nThreads()-1; ++i) {
            executables.push_back(executor.pack(worker, it, it + grainSize));
            it += grainSize;
        }
        executables.push_back(executor.pack(worker, it, cix.size()));
        executor.execute_and_wait(std::move(executables));
    }
}

void ContiguousCLLNeighborList::update(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if(!_is_set_up) {
        set_up(node.subnode("setUp"));
    } else {
        ccll.update(node.subnode("update CLL"));
        fill_verlet_list(node.subnode("fill verlet list"));
    }
}

void ContiguousCLLNeighborList::clear(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (_max_cutoff > 0) {
        ccll.clear();
        for (auto &neighbors : _data.get().neighbors) {
            neighbors.clear();
        }
    }
}

void ContiguousCLLNeighborList::updateData(model::CPUParticleData::update_t &&update) {
    _data.get().update(std::forward<model::CPUParticleData::update_t>(update));
}


}
}
}
}