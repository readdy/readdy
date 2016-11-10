/**
 * << detailed description >>
 *
 * @file NeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 08.09.16
 * @todo hilbert curve sorting of cells, grouping of particles into cells, replace maps by vector w/ padding
 */

#include <readdy/kernel/cpu/model/NeighborList.h>
#include <sstream>
#include <readdy/kernel/cpu/util/hilbert.h>

template<typename T, typename Predicate>
typename std::vector<T>::iterator insert_sorted(std::vector<T> &vec, T const &item, Predicate pred) {
    return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, pred), item);
}

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

template<typename T, typename Dims>
T get_contiguous_index(T i, T j, T k, Dims I, Dims J) {
    return k + j * J + i * I * J;
};

auto cells_insert_predicate = [](const NeighborList::Cell& c1, const NeighborList::Cell& c2) -> bool {
    return c1.hilbert_index < c2.hilbert_index;
};

auto cells_search_predicate = [](const NeighborList::Cell& c1, const NeighborList::hilbert_index_t idx) -> bool {
    return c1.hilbert_index < idx;
};

const static std::vector<NeighborList::neighbor_t> no_neighbors{};

NeighborList::NeighborList(const readdy::model::KernelContext *const context, util::Config const *const config,
                           skin_size_t skin)
        : ctx(context), config(config), cells(std::vector<Cell>()), simBoxSize(ctx->getBoxSize()),
          maps(config->nThreads), skin_size(skin) {}

void NeighborList::setupCells() {
    if (cells.empty()) {
        double maxCutoff = 0;
        for (auto &&e : ctx->getAllOrder2RegisteredPotentialTypes()) {
            for (auto &&p : ctx->getOrder2Potentials(std::get<0>(e), std::get<1>(e))) {
                maxCutoff = maxCutoff < p->getCutoffRadius() ? p->getCutoffRadius() : maxCutoff;
            }
        }
        for (auto &&e : ctx->getAllOrder2Reactions()) {
            maxCutoff = maxCutoff < e->getEductDistance() ? e->getEductDistance() : maxCutoff;
        }
        NeighborList::maxCutoff = maxCutoff + skin_size;
        if (maxCutoff > 0) {
            // todo sort cells by hilbert curve
            const auto desiredCellWidth = .5 * maxCutoff;

            for (unsigned short i = 0; i < 3; ++i) {
                nCells[i] = static_cast<cell_index>(floor(simBoxSize[i] / desiredCellWidth));
                if (nCells[i] == 0) nCells[i] = 1;
                cellSize[i] = simBoxSize[i] / static_cast<double>(nCells[i]);
            }
            /**
             * const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
             * const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
             * const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
             *
             * const auto &periodic = ctx->getPeriodicBoundary();
             * if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
             * else if (i < 0 || i >= nCells[0]) return nullptr;
             * if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
             * else if (j < 0 || j >= nCells[1]) return nullptr;
             * if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
             * else if (k < 0 || k >= nCells[2]) return nullptr;
             * return &cells.at(k + j * nCells[2] + i * nCells[2] * nCells[1]);
             */

            hilbertIndexMapping.resize(nCells[0] * nCells[1] * nCells[2]);
            auto it_mapping = hilbertIndexMapping.begin();
            for (cell_index i = 0; i < nCells[0]; ++i) {
                for (cell_index j = 0; j < nCells[1]; ++j) {
                    for (cell_index k = 0; k < nCells[2]; ++k) {
                        auto it = insert_sorted(cells, {i, j, k, nCells}, cells_insert_predicate);
                        *it_mapping = it->hilbert_index;
                        ++it_mapping;
                    }
                }
            }

            for (cell_index i = 0; i < nCells[0]; ++i) {
                for (cell_index j = 0; j < nCells[1]; ++j) {
                    for (cell_index k = 0; k < nCells[2]; ++k) {
                        setupNeighboringCells(i, j, k);
                    }
                }
            }
        } else {
            nCells = {{1, 1, 1}};
            cells.push_back({0, 0, 0, nCells});
            hilbertIndexMapping.push_back(cells.back().hilbert_index);
            setupNeighboringCells(0, 0, 0);
        }
    }
    if (cells.size() > 1) {
        // todo group particles by cells
    }
}

void
NeighborList::setupNeighboringCells(const signed_cell_index i, const signed_cell_index j, const signed_cell_index k) {
    auto me = getCell(i, j, k);
    for (signed_cell_index _i = -2; _i < 3; ++_i) {
        for (signed_cell_index _j = -2; _j < 3; ++_j) {
            for (signed_cell_index _k = -2; _k < 3; ++_k) {
                // don't add me as neighbor to myself
                if (!(_i == 0 && _j == 0 && _k == 0)) {
                    me->addNeighbor(getCell(i + _i, j + _j, k + _k));
                }
            }
        }
    }
}

void NeighborList::clear() {
    cells.clear();
    maps.clear();
}

void NeighborList::fillCells(data_t &data) {
    if (maxCutoff > 0) {

        for (auto &cell : cells) {
            cell.particleIndices.clear();
        }
        std::for_each(maps.begin(), maps.end(), [](container_t &map) { map.clear(); });
        maps.resize(config->nThreads);

        auto it = data.entries.begin();
        while (it != data.entries.end()) {
            if (!it->is_deactivated()) {
                auto cell = getCell(it->pos);
                if (cell) {
                    auto ptr = &(*it);
                    getPairs(cell).emplace(ptr, std::vector<neighbor_t>());
                    cell->particleIndices.push_back(ptr);
                }
            }
            ++it;
        }

        {
            using cell_it_t = decltype(cells.begin());
            const auto size = cells.size();
            const std::size_t grainSize = size / config->nThreads;
            const auto cutoffSquared = maxCutoff * maxCutoff;
            auto worker = [&data, cutoffSquared, this](const cell_it_t begin, const cell_it_t end, container_t &map) {
                const auto &d2 = ctx->getDistSquaredFun();
                for (auto _b = begin; _b != end; ++_b) {
                    auto &cell = *_b;
                    for (const auto &pI : cell.particleIndices) {
                        const auto pIIt = map.find(pI);
                        if (pIIt != map.end()) {
                            auto &pI_vec = pIIt->second;
                            for (const auto &pJ : cell.particleIndices) {
                                if (pI != pJ) {
                                    const auto distSquared = d2(pI->pos, pJ->pos);
                                    if (distSquared < cutoffSquared) {
                                        Neighbor neighbor{pJ, distSquared};
                                        pI_vec.push_back(std::move(neighbor));
                                    }
                                }
                            }
                            for (auto &neighboringCell : cell.neighbors) {
                                for (const auto pJ : neighboringCell->particleIndices) {
                                    const auto distSquared = d2(pI->pos, pJ->pos);
                                    if (distSquared < cutoffSquared) {
                                        Neighbor neighbor{pJ, distSquared};
                                        pI_vec.push_back(std::move(neighbor));
                                    }
                                }
                            }

                        }
                    }
                }
            };

            std::vector<util::scoped_thread> threads;
            threads.reserve(config->nThreads);

            auto it_cells = cells.begin();
            auto it_maps = maps.begin();
            while (it_maps != maps.end() - 1) {
                threads.push_back(
                        util::scoped_thread(std::thread(worker, it_cells, it_cells + grainSize, std::ref(*it_maps))));
                it_cells += grainSize;
                ++it_maps;
            }
            threads.push_back(util::scoped_thread(std::thread(worker, it_cells, cells.end(), std::ref(*it_maps))));
        }
    }
}

void NeighborList::create(data_t &data) {
    simBoxSize = ctx->getBoxSize();
    setupCells();
    fillCells(data);
}

NeighborList::Cell *NeighborList::getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k) {
    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
    else if (i < 0 || i >= nCells[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
    else if (j < 0 || j >= nCells[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
    else if (k < 0 || k >= nCells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, nCells[1], nCells[2]);
    if(cix >= hilbertIndexMapping.size()) {
        log::console()->error("shit squared2!");
    }
    const auto hix = hilbertIndexMapping.at(cix);
    if(hix > cells.back().hilbert_index) {
        log::console()->error("shit");
    }
    return &*std::lower_bound(cells.begin(), cells.end(), hix, cells_search_predicate);
}

void NeighborList::remove(const particle_index idx) {
    auto cell = getCell(idx->pos);
    if (cell != nullptr) {
        std::set<container_t *> affectedMaps;
        for (auto neighborCell : getCell(idx->pos)->neighbors) {
            affectedMaps.insert(&getPairs(neighborCell));
        }
        try {
            for (auto &neighbor : neighbors(idx)) {
                auto &&neighborCell = getCell(neighbor.idx->pos);
                for (auto map : affectedMaps) {
                    const auto it = map->find(neighbor.idx);
                    if (it != map->end()) {
                        std::remove_if(it->second.begin(), it->second.end(),
                                       [idx](const neighbor_t &n) {
                                           return n.idx == idx;
                                       });
                    }

                }
            }
            getPairs(cell).erase(idx);
        } catch (const std::out_of_range &) {
            log::console()->error("tried to remove particle with id {} but it was not in the neighbor list", idx->id);
        }
    }
}

void NeighborList::insert(const data_t &data, const particle_index idx) {
    const auto d2 = ctx->getDistSquaredFun();
    const auto pos = idx->pos;
    const auto cutoffSquared = maxCutoff * maxCutoff;
    auto cell = getCell(pos);
    if (cell) {
        auto &map = getPairs(cell);
        cell->particleIndices.push_back(idx);
        auto emplace_ret = map.emplace(idx, std::vector<neighbor_t>());

        for (const auto pJ : cell->particleIndices) {
            if (idx != pJ) {
                const auto distSquared = d2(pos, pJ->pos);
                if (distSquared < cutoffSquared) {
                    (*(emplace_ret.first)).second.push_back({pJ, distSquared});
                    map[pJ].push_back({idx, distSquared});
                }
            }
        }
        for (auto &neighboringCell : cell->neighbors) {
            auto &neighboring_map = getPairs(neighboringCell);
            for (const auto &pJ : neighboringCell->particleIndices) {
                const auto distSquared = d2(pos, pJ->pos);
                if (distSquared < cutoffSquared) {
                    (*(emplace_ret.first)).second.push_back({pJ, distSquared});
                    neighboring_map[pJ].push_back({idx, distSquared});
                }
            }
        }
    } else {
        //log::console()->error("could not assign particle (index={}) to any cell!", data.getEntryIndex(idx));
    }
}

NeighborList::Cell *NeighborList::getCell(const readdy::model::Particle::pos_type &pos) {
    const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    return getCell(i, j, k);
}

void NeighborList::updateData(NeighborList::data_t &data, ParticleData::update_t update) {
    if (maxCutoff > 0) {
        for (const auto &p : std::get<1>(update)) {
            remove(p);
        }
    }

    auto newEntries = data.update(std::move(update));
    if (maxCutoff > 0) {
        for (const auto &p : newEntries) {
            insert(data, p);
        }
    }
}

const std::vector<NeighborList::neighbor_t> &NeighborList::neighbors(NeighborList::particle_index const entry) const {
    if (maxCutoff > 0) {
        const auto cell = getCell(entry->pos);
        if (cell != nullptr) {
            return getPairs(cell).at(entry);
        } else {
            std::stringstream stream;
            stream << "the particle position " << entry->pos << " was in no neighbor list cell" << std::endl;
            throw std::out_of_range(stream.str());
        }
    }
    return no_neighbors;
}

const NeighborList::Cell *const NeighborList::getCell(const readdy::model::Particle::pos_type &pos) const {
    const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    return getCell(i, j, k);
}

const NeighborList::Cell *const NeighborList::getCell(NeighborList::signed_cell_index i,
                                                      NeighborList::signed_cell_index j,
                                                      NeighborList::signed_cell_index k) const {

    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
    else if (i < 0 || i >= nCells[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
    else if (j < 0 || j >= nCells[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
    else if (k < 0 || k >= nCells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, nCells[1], nCells[2]);
    if(cix >= hilbertIndexMapping.size()) {
        log::console()->error("shit squared!");
    }
    const auto hix = hilbertIndexMapping.at(cix);
    if(hix > cells.back().hilbert_index) {
        log::console()->error("shit");
    }
    return &*std::lower_bound(cells.begin(), cells.end(), hix, cells_search_predicate);
}

NeighborList::container_t &NeighborList::getPairs(const NeighborList::Cell *const cell) {
    return maps.at(getMapsIndex(cell));
}

const NeighborList::container_t &NeighborList::getPairs(const NeighborList::Cell *const cell) const {
    return maps.at(getMapsIndex(cell));
}

util::Config::n_threads_t NeighborList::getMapsIndex(const NeighborList::Cell *const cell) const {
    return static_cast<util::Config::n_threads_t>(floor(cell->contiguous_index / floor(cells.size() / config->nThreads)));
}

const std::vector<NeighborList::neighbor_t> &NeighborList::find_neighbors(particle_index const entry) const {
    if (maxCutoff > 0) {
        const auto cell = getCell(entry->pos);
        if (cell != nullptr) {
            const auto &map = getPairs(cell);
            const auto it = map.find(entry);
            if (it != map.end()) {
                return it->second;
            }
        }
    }
    return no_neighbors;
}

void NeighborList::displace(NeighborList::data_iter_t iter, const readdy::model::Vec3 &vec) {
    iter->pos += vec;
    // todo record displacement
}

NeighborList::~NeighborList() = default;

Neighbor::Neighbor(const index_t idx, const double d2) : idx(idx), d2(d2) {}

Neighbor::Neighbor(Neighbor &&rhs) : idx(rhs.idx), d2(std::move(rhs.d2)) {
}

Neighbor &Neighbor::operator=(Neighbor &&rhs) {
    idx = rhs.idx;
    d2 = std::move(rhs.d2);
    return *this;
}

void NeighborList::Cell::addNeighbor(NeighborList::Cell *cell) {
    if (cell && cell->hilbert_index != hilbert_index
        && (enoughCells || std::find(neighbors.begin(), neighbors.end(), cell) == neighbors.end())) {
        neighbors.push_back(cell);
    }
}

NeighborList::Cell::Cell(cell_index i, cell_index j, cell_index k,
                         const std::array<NeighborList::cell_index, 3> &nCells)
        : hilbert_index(getHilbertIndex(i, j, k)),
          contiguous_index(get_contiguous_index(i, j, k, nCells[1], nCells[2])),
          enoughCells(nCells[0] >= 5 && nCells[1] >= 5 && nCells[2] >= 5) {
}

bool operator==(const NeighborList::Cell &lhs, const NeighborList::Cell &rhs) {
    return lhs.contiguous_index == rhs.contiguous_index;
}

bool operator!=(const NeighborList::Cell &lhs, const NeighborList::Cell &rhs) {
    return !(lhs == rhs);
}

NeighborList::hilbert_index_t getHilbertIndex(NeighborList::cell_index i,
                                              NeighborList::cell_index j,
                                              NeighborList::cell_index k) {
    bitmask_t coords[3] {i,j,k};
    return static_cast<unsigned int>(hilbert_c2i(3, CHAR_BIT, coords));
}

}
}
}
}
