/**
 * << detailed description >>
 *
 * @file NeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 08.09.16
 * @todo hilbert curve sorting of cells, grouping of particles into cells, replace maps by vector w/ padding
 */

#include <set>
#include <sstream>

#include <readdy/kernel/cpu/model/NeighborList.h>
#include <readdy/kernel/cpu/util/hilbert.h>
#include <readdy/common/numeric.h>

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

const static std::vector<NeighborList::neighbor_t> no_neighbors{};

NeighborList::NeighborList(const ctx_t *const context, data_t &data, util::Config const *const config,
                           skin_size_t skin)
        : ctx(context), config(config), cells(std::vector<Cell>()),
          simBoxSize(ctx->getBoxSize()), skin_size(skin), data(data) {}

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

            hilbertIndexMapping.reserve(nCells[0] * nCells[1] * nCells[2]);
            for (cell_index i = 0; i < nCells[0]; ++i) {
                for (cell_index j = 0; j < nCells[1]; ++j) {
                    for (cell_index k = 0; k < nCells[2]; ++k) {
                        cells.push_back({i, j, k, nCells});
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
    data.neighbors.clear();
}

void NeighborList::fillCells() {
    if (maxCutoff > 0) {

        for (auto &cell : cells) {
            cell.particleIndices.clear();
        }
        data.neighbors.clear();
        data.neighbors.resize(data.size() + data.getNDeactivated());

        data_t::index_t idx = 0;
        for(const auto& entry : data) {
            if (!entry.is_deactivated()) {
                auto cell = getCell(entry.pos);
                if (cell) {
                    cell->particleIndices.push_back(idx);
                }
            }
            ++idx;
        }

        {
            using cell_it_t = decltype(cells.begin());
            const auto size = cells.size();
            const std::size_t grainSize = size / config->nThreads;
            const auto cutoffSquared = maxCutoff * maxCutoff;
            auto worker = [cutoffSquared, this](const cell_it_t begin, const cell_it_t end) {
                const auto &d2 = ctx->getDistSquaredFun();
                for (cell_it_t _b = begin; _b != end; ++_b) {
                    auto &cell = *_b;
                    for (const auto &pI : cell.particleIndices) {
                        const auto& particle_i = data.entries.at(pI);
                        auto& neighbors_i = data.neighbors.at(pI);
                        for (const auto &pJ : cell.particleIndices) {
                            if (pI != pJ) {
                                const auto distSquared = d2(particle_i.pos, data.pos(pJ));
                                if (distSquared < cutoffSquared) {
                                    neighbors_i.push_back({pJ, distSquared});
                                }
                            }
                        }
                        for (auto &neighboringCell : cell.neighbors) {
                            for (const auto pJ : neighboringCell->particleIndices) {
                                const auto distSquared = d2(particle_i.pos, data.pos(pJ));
                                if (distSquared < cutoffSquared) {
                                    neighbors_i.push_back({pJ, distSquared});
                                }
                            }
                        }

                    }
                }
            };

            std::vector<util::scoped_thread> threads;
            threads.reserve(config->nThreads);

            auto it_cells = cells.begin();
            for(int i = 0; i < config->nThreads-1; ++i) {
                threads.push_back(
                        util::scoped_thread(std::thread(worker, it_cells, it_cells + grainSize)));
                it_cells += grainSize;
            }
            threads.push_back(util::scoped_thread(std::thread(worker, it_cells, cells.end())));
        }
    }
}

void NeighborList::create() {
    simBoxSize = ctx->getBoxSize();
    data.setFixPosFun(ctx->getFixPositionFun());
    setupCells();
    fillCells();
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
    return &cells.at(static_cast<cell_index>(cix));
}

void NeighborList::remove(const particle_index idx) {
    auto cell = getCell(data.pos(idx));
    if (cell != nullptr) {
        auto remove_predicate = [idx] (const ParticleData::Neighbor& n){
            return n.idx == idx;
        };
        try {
            for (auto &neighbor : neighbors(idx)) {
                auto& neighbors_2nd = data.neighbors.at(neighbor.idx);
                neighbors_2nd.erase(std::find_if(neighbors_2nd.begin(), neighbors_2nd.end(), remove_predicate));
            }
            auto find_it = std::find(cell->particleIndices.begin(), cell->particleIndices.end(), idx);
            if(find_it != cell->particleIndices.end()) {
                cell->particleIndices.erase(find_it);
            }
        } catch (const std::out_of_range &) {
            log::console()->error("tried to remove particle with id {} but it was not in the neighbor list", idx);
        }
    }
}

void NeighborList::insert(const particle_index idx) {
    const auto& d2 = ctx->getDistSquaredFun();
    const auto& pos = data.pos(idx);
    const auto cutoffSquared = maxCutoff * maxCutoff;
    auto cell = getCell(pos);
    if (cell) {
        cell->particleIndices.push_back(idx);
        auto& myNeighbors = data.neighbors.at(idx);
        for (const auto pJ : cell->particleIndices) {
            if (idx != pJ) {
                const auto distSquared = d2(pos, data.pos(pJ));
                if (distSquared < cutoffSquared) {
                    myNeighbors.push_back({pJ, distSquared});
                    data.neighbors.at(pJ).push_back({idx, distSquared});
                }
            }
        }
        for (auto &neighboringCell : cell->neighbors) {
            for (const auto &pJ : neighboringCell->particleIndices) {
                const auto distSquared = d2(pos, data.pos(pJ));
                if (distSquared < cutoffSquared) {
                    myNeighbors.push_back({pJ, distSquared});
                    data.neighbors.at(pJ).push_back({idx, distSquared});
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

void NeighborList::updateData(ParticleData::update_t update) {
    if (maxCutoff > 0) {
        for (const auto &p : std::get<1>(update)) {
            remove(p);
        }
    }

    auto newEntries = data.update(std::move(update));
    if (maxCutoff > 0) {
        for (const auto &p : newEntries) {
            insert(p);
        }
    }
}

const std::vector<NeighborList::neighbor_t> &NeighborList::neighbors(NeighborList::particle_index const entry) const {
    if (maxCutoff > 0) {
        return data.neighbors.at(entry);
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
    return &cells.at(cix);
}

util::Config::n_threads_t NeighborList::getMapsIndex(const NeighborList::Cell *const cell) const {
    return static_cast<util::Config::n_threads_t>(floor(cell->contiguous_index / floor(cells.size() / config->nThreads)));
}

const std::vector<NeighborList::neighbor_t> &NeighborList::find_neighbors(particle_index const entry) const {
    if (maxCutoff > 0 && entry < data.neighbors.size()) {
        return data.neighbors.at(entry);
    }
    return no_neighbors;
}

void NeighborList::displace(NeighborList::data_iter_t iter, const readdy::model::Vec3 &vec) {
    data.displace(*iter, vec);
}

NeighborList::~NeighborList() = default;

void NeighborList::Cell::addNeighbor(NeighborList::Cell *cell) {
    if (cell && cell->contiguous_index != contiguous_index
        && (enoughCells || std::find(neighbors.begin(), neighbors.end(), cell) == neighbors.end())) {
        neighbors.push_back(cell);
    }
}

NeighborList::Cell::Cell(cell_index i, cell_index j, cell_index k,
                         const std::array<NeighborList::cell_index, 3> &nCells)
        :  contiguous_index(get_contiguous_index(i, j, k, nCells[1], nCells[2])),
          enoughCells(nCells[0] >= 5 && nCells[1] >= 5 && nCells[2] >= 5) {
}

bool operator==(const NeighborList::Cell &lhs, const NeighborList::Cell &rhs) {
    return lhs.contiguous_index == rhs.contiguous_index;
}

bool operator!=(const NeighborList::Cell &lhs, const NeighborList::Cell &rhs) {
    return !(lhs == rhs);
}

NeighborList::hilbert_index_t NeighborList::getHilbertIndex(std::size_t i, std::size_t j, std::size_t k) const {
    bitmask_t coords[3] {i,j,k};
    return static_cast<unsigned int>(hilbert_c2i(3, CHAR_BIT, coords));
}


}
}
}
}
