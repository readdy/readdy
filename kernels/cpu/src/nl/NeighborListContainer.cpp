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
 * @file NeighborListContainer.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */


#include <readdy/kernel/cpu/nl/NeighborListContainer.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {


NeighborListContainer::NeighborListContainer(NLContainerConfig config)
        : _indexer(), _config(config), _elements() {}

NeighborListContainer::const_iterator NeighborListContainer::begin() const {
    return _elements.begin();
}

typename NeighborListContainer::const_iterator NeighborListContainer::end() const {
    return _elements.end();
}

void NeighborListContainer::clear() {
    _elements.clear();
    _indexer.clear();
}

void ContiguousCLLNeighborListContainer::update(scalar cutoffSquared, const util::PerformanceNode &perf) {
    auto t = perf.timeit();
    if(_elements.size() != _config.threads.get().nThreads()) {
        _elements.clear();
        _elements.resize(_config.threads.get().nThreads());
    } else {
        std::for_each(_elements.begin(), _elements.end(), [](local_index_vector &vec) {
            vec.resize(0);
        });
    }

    const auto &cll = _cll;

    auto cix = cll.cellIndex();

    const auto grainSize = cix.size() / _config.threads.get().nThreads();

    const auto &data = _config.data.get();
    const auto &context = _config.context.get();
    const auto &d2 = context.distSquaredFun();

    auto worker = [this, cix, &data, &cll, &d2, cutoffSquared](std::size_t tid, std::size_t begin, std::size_t end) {
        auto &neighbors = _elements.at(tid);

        for (std::size_t cellIndex = begin; cellIndex < end; ++cellIndex) {

            for (auto itParticles = cll.particlesBegin(cellIndex);
                 itParticles != cll.particlesEnd(cellIndex); ++itParticles) {

                auto pidx = *itParticles;
                const auto &entry = data.entry_at(pidx);

                if(!entry.deactivated) {

                    neighbors.push_back(pidx);

                    auto n_neighbors_index = neighbors.size();
                    neighbors.push_back(0_z);
                    auto n_neighbors = 0_z;

                    for (auto itPP = cll.particlesBegin(cellIndex); itPP != cll.particlesEnd(cellIndex); ++itPP) {
                        auto pparticle = *itPP;
                        if (pidx != pparticle) {
                            const auto &pp = data.entry_at(pparticle);
                            if (!pp.deactivated) {
                                if (d2(entry.pos, pp.pos) < cutoffSquared) {
                                    neighbors.push_back(pparticle);
                                    ++n_neighbors;
                                }
                            }
                        }
                    }

                    for (auto itNeighborCell = cll.neighborsBegin(cellIndex);
                         itNeighborCell != cll.neighborsEnd(cellIndex); ++itNeighborCell) {
                        for (auto itNeighborParticle = cll.particlesBegin(*itNeighborCell);
                             itNeighborParticle != cll.particlesEnd(*itNeighborCell); ++itNeighborParticle) {
                            const auto &neighbor = data.entry_at(*itNeighborParticle);
                            if (!neighbor.deactivated) {
                                if (d2(entry.pos, neighbor.pos) < cutoffSquared) {
                                    neighbors.push_back(*itNeighborParticle);
                                    ++n_neighbors;
                                }
                            }
                        }
                    }

                    if (n_neighbors > 0) neighbors.at(n_neighbors_index) = n_neighbors;
                }
            }
        }
    };
    const auto &executor = *_config.threads.get().executor();
    std::vector<std::function<void(std::size_t)>> executables;
    executables.reserve(_config.threads.get().nThreads());
    auto it = 0_z;
    for (int i = 0; i < _config.threads.get().nThreads() - 1; ++i) {
        executables.push_back(executor.pack(worker, it, it + grainSize));
        it += grainSize;
    }
    executables.push_back(executor.pack(worker, it, cix.size()));
    executor.execute_and_wait(std::move(executables));
}


NLContainerConfig::NLContainerConfig(const model::KernelContext &context,
                                     const NLContainerConfig::thread_config_type &config,
                                     const NLContainerConfig::data_container_type &data)
        : context(context), threads(config), data(data) {
}

ContiguousCLLNeighborListContainer::ContiguousCLLNeighborListContainer(NLContainerConfig config,
                                                                       const ContiguousCellLinkedList &cll)
        : NeighborListContainer(config), _cll(cll) {}

DynamicCLLNeighborListContainer::DynamicCLLNeighborListContainer(NLContainerConfig config,
                                                                 const DynamicCellLinkedList &cll)
        : NeighborListContainer(config), _cll(cll) {}

void DynamicCLLNeighborListContainer::update(scalar cutoffSquared, const util::PerformanceNode &perf) {
    auto t = perf.timeit();
    if(_elements.size() != _config.threads.get().nThreads()) {
        _elements.clear();
        _elements.resize(_config.threads.get().nThreads());
    } else {
        std::for_each(_elements.begin(), _elements.end(), [](local_index_vector &vec) {
            vec.resize(0);
        });
    }

    auto cix = _cll.cellIndex();
    const auto &d2 = _config.context.get().distSquaredFun();


    const auto grainSize = cix.size() / _config.threads.get().nThreads();
    const auto &data = _config.data.get();
    const auto &dcll = _cll;
    auto worker = [this, cix, &dcll, &d2, &data, cutoffSquared](std::size_t tid, std::size_t begin, std::size_t end) {

        auto &neighbors = _elements.at(tid);
        for (std::size_t cellIndex = begin; cellIndex < end; ++cellIndex) {
            for (auto itParticles = dcll.particlesBegin(cellIndex);
                 itParticles != dcll.particlesEnd(cellIndex); ++itParticles) {
                auto pidx = *itParticles;
                auto &entry = data.entry_at(pidx);
                if(!entry.deactivated) {
                    neighbors.push_back(pidx);
                    auto n_neighbors_index = neighbors.size();
                    neighbors.push_back(0_z);
                    auto n_neighbors = 0_z;


                    for (auto itPP = dcll.particlesBegin(cellIndex); itPP != dcll.particlesEnd(cellIndex); ++itPP) {
                        auto pparticle = *itPP;
                        if (pidx != pparticle) {
                            const auto &pp = data.entry_at(pparticle);
                            if (!pp.deactivated) {
                                const auto distSquared = d2(entry.pos, pp.pos);
                                if (distSquared < cutoffSquared) {
                                    neighbors.push_back(pparticle);
                                    ++n_neighbors;
                                }
                            }
                        }
                    }

                    for (auto itNeighborCell = dcll.neighborsBegin(cellIndex);
                         itNeighborCell != dcll.neighborsEnd(cellIndex); ++itNeighborCell) {
                        for (auto itNeighborParticle = dcll.particlesBegin(*itNeighborCell);
                             itNeighborParticle != dcll.particlesEnd(*itNeighborCell); ++itNeighborParticle) {
                            const auto &neighbor = data.entry_at(*itNeighborParticle);
                            if (!neighbor.deactivated) {
                                const auto distSquared = d2(entry.pos, neighbor.pos);
                                if (distSquared < cutoffSquared) {
                                    neighbors.push_back(*itNeighborParticle);
                                    ++n_neighbors;
                                }
                            }
                        }
                    }
                    if (n_neighbors > 0) neighbors.at(n_neighbors_index) = n_neighbors;
                }
            }
        }
    };
    const auto &executor = *_config.threads.get().executor();
    std::vector<std::function<void(std::size_t)>> executables;
    executables.reserve(_config.threads.get().nThreads());
    auto it = 0_z;
    for (int i = 0; i < _config.threads.get().nThreads() - 1; ++i) {
        executables.push_back(executor.pack(worker, it, it + grainSize));
        it += grainSize;
    }
    executables.push_back(executor.pack(worker, it, cix.size()));
    executor.execute_and_wait(std::move(executables));
}


void CompactCLLNeighborListContainer::update(scalar cutoffSquared, const util::PerformanceNode &perf) {
    if(serialUpdate) {
        updateSerial(cutoffSquared);
    } else {
        updateParallel(cutoffSquared, perf);
    }
}

CompactCLLNeighborListContainer::CompactCLLNeighborListContainer(NLContainerConfig config,
                                                                 const CompactCellLinkedList &cll)
        : NeighborListContainer(config), _cll(cll) {}

void CompactCLLNeighborListContainer::updateSerial(scalar cutoffSquared) {
    if(_elements.size() != _config.threads.get().nThreads()) {
        _elements.clear();
        _elements.resize(_config.threads.get().nThreads());
    }

    auto cix = _cll.cellIndex();

    const auto &d2 = _config.context.get().distSquaredFun();
    const auto &data = _config.data.get();
    const auto &ccll = _cll;

    const auto &head = ccll.head();
    const auto &list = ccll.list();

    _elements.resize(1);
    auto &neighbors = _elements.at(0);
    neighbors.clear();

    for (std::size_t cellIndex = 0; cellIndex < cix.size(); ++cellIndex) {

        auto pptr = (*head.at(cellIndex)).load();
        while (pptr != 0) {
            auto pidx = pptr - 1;
            const auto &entry = data.entry_at(pidx);
            if(!entry.deactivated) {
                neighbors.push_back(pidx);
                auto n_neighbors_index = neighbors.size();
                neighbors.push_back(0_z);
                auto n_neighbors = 0_z;

                {
                    auto ppptr = (*head.at(cellIndex)).load();
                    while (ppptr != 0) {
                        auto ppidx = ppptr - 1;
                        if (ppidx != pidx) {
                            const auto &pp = data.entry_at(ppidx);
                            if (!pp.deactivated) {
                                const auto distSquared = d2(entry.pos, pp.pos);
                                if (distSquared < cutoffSquared) {
                                    neighbors.push_back(ppidx);
                                    ++n_neighbors;
                                }
                            }
                        }
                        ppptr = list.at(ppptr);
                    }
                }

                for (auto itNeighborCell = ccll.neighborsBegin(cellIndex);
                     itNeighborCell != ccll.neighborsEnd(cellIndex); ++itNeighborCell) {

                    auto nptr = (*head.at(*itNeighborCell)).load();
                    while (nptr != 0) {
                        auto nidx = nptr - 1;

                        const auto &neighbor = data.entry_at(nidx);
                        if (!neighbor.deactivated) {
                            const auto distSquared = d2(entry.pos, neighbor.pos);
                            if (distSquared < cutoffSquared) {
                                neighbors.push_back(nidx);
                                ++n_neighbors;
                            }
                        }

                        nptr = list.at(nptr);
                    }
                }
                if (n_neighbors > 0) neighbors.at(n_neighbors_index) = n_neighbors;
            }
            pptr = list.at(pptr);
        }
    }
}

void CompactCLLNeighborListContainer::updateParallel(scalar cutoffSquared, const util::PerformanceNode &perf) {
    auto t = perf.timeit();
    if(_elements.size() != _config.threads.get().nThreads()) {
        _elements.clear();
        _elements.resize(_config.threads.get().nThreads());
    } else {
        std::for_each(_elements.begin(), _elements.end(), [](local_index_vector &vec) {
            vec.resize(0);
        });
    }

    auto cix = _cll.cellIndex();

    const auto &d2 = _config.context.get().distSquaredFun();
    const auto &data = _config.data.get();
    const auto &ccll = _cll;

    const auto grainSize = cix.size() / _config.threads.get().nThreads();
    auto worker = [this, cix, &d2, &data, &ccll, cutoffSquared](std::size_t tid, std::size_t begin, std::size_t end) {
        const auto &head = ccll.head();
        const auto &list = ccll.list();

        auto &neighbors = _elements.at(tid);
        for (std::size_t cellIndex = begin; cellIndex < end; ++cellIndex) {

            auto pptr = (*head.at(cellIndex)).load();
            while (pptr != 0) {
                auto pidx = pptr - 1;
                const auto &entry = data.entry_at(pidx);
                if(!entry.deactivated) {
                    neighbors.push_back(pidx);
                    auto n_neighbors_index = neighbors.size();
                    neighbors.push_back(0_z);
                    auto n_neighbors = 0_z;

                    {
                        auto ppptr = (*head.at(cellIndex)).load();
                        while (ppptr != 0) {
                            auto ppidx = ppptr - 1;
                            if (ppidx != pidx) {
                                const auto &pp = data.entry_at(ppidx);
                                if (!pp.deactivated) {
                                    const auto distSquared = d2(entry.pos, pp.pos);
                                    if (distSquared < cutoffSquared) {
                                        neighbors.push_back(ppidx);
                                        ++n_neighbors;
                                    }
                                }
                            }
                            ppptr = list.at(ppptr);
                        }
                    }

                    for (auto itNeighborCell = ccll.neighborsBegin(cellIndex);
                         itNeighborCell != ccll.neighborsEnd(cellIndex); ++itNeighborCell) {

                        auto nptr = (*head.at(*itNeighborCell)).load();
                        while (nptr != 0) {
                            auto nidx = nptr - 1;

                            const auto &neighbor = data.entry_at(nidx);
                            if (!neighbor.deactivated) {
                                const auto distSquared = d2(entry.pos, neighbor.pos);
                                if (distSquared < cutoffSquared) {
                                    neighbors.push_back(nidx);
                                    ++n_neighbors;
                                }
                            }

                            nptr = list.at(nptr);
                        }
                    }
                    if(n_neighbors > 0) neighbors.at(n_neighbors_index) = n_neighbors;
                }
                pptr = list.at(pptr);
            }
        }
    };
    const auto &executor = *_config.threads.get().executor();
    std::vector<std::function<void(std::size_t)>> executables;
    executables.reserve(_config.threads.get().nThreads());
    auto it = 0_z;
    for (int i = 0; i < _config.threads.get().nThreads() - 1; ++i) {
        executables.push_back(executor.pack(worker, it, it + grainSize));
        it += grainSize;
    }
    executables.push_back(executor.pack(worker, it, cix.size()));
    executor.execute_and_wait(std::move(executables));
}


}
}
}
}