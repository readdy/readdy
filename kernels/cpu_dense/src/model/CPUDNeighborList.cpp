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
 * @date 23.11.16
 */

#include <readdy/common/numeric.h>
#include <readdy/common/thread/scoped_thread.h>
#include <readdy/common/thread/barrier.h>
#include "readdy/kernel/cpu_dense/model/CPUDNeighborList.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace model {


CPUDNeighborList::Neighbor::Neighbor(const index_t idx, const double d2) : idx(idx), d2(d2) {
}

namespace thd = readdy::util::thread;

template<typename T, typename Dims>
T get_contiguous_index(T i, T j, T k, Dims I, Dims J) {
    return k + j * J + i * I * J;
};

const static std::vector<CPUDNeighborList::neighbor_t> no_neighbors{};

CPUDNeighborList::CPUDNeighborList(const ctx_t *const context, data_t &data, readdy::util::thread::Config const *const config)
        : ctx(context), config(config), cells(std::vector<Cell>()), simBoxSize(ctx->getBoxSize()), data(data) {}

void CPUDNeighborList::setupCells() {
    if (cells.empty()) {
        double maxCutoff = 0;
        for(const auto& entry : ctx->potentials().potentials_order2()) {
            for(const auto& potential : entry.second) {
                maxCutoff = maxCutoff < potential->getCutoffRadius() ? potential->getCutoffRadius() : maxCutoff;
            }
        }
        for (auto &&e : ctx->reactions().order2_flat()) {
            maxCutoff = maxCutoff < e->getEductDistance() ? e->getEductDistance() : maxCutoff;
        }
        CPUDNeighborList::maxCutoff = maxCutoff;

        if (maxCutoff > 0) {
            const auto desiredCellWidth = .5 * maxCutoff;

            for (unsigned short i = 0; i < 3; ++i) {
                nCells[i] = static_cast<cell_index>(floor(simBoxSize[i] / desiredCellWidth));
                if (nCells[i] == 0) nCells[i] = 1;
                cellSize[i] = simBoxSize[i] / static_cast<double>(nCells[i]);
            }
            log::debug("resulting cell size = {}", cellSize);

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
}

void
CPUDNeighborList::setupNeighboringCells(const signed_cell_index i, const signed_cell_index j, const signed_cell_index k) {
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

void CPUDNeighborList::clear() {
    cells.clear();
    for (auto &list : neighbor_list) {
        list.clear();
    }
}

void CPUDNeighborList::fillCells() {
    if (maxCutoff > 0) {
        const auto c2 = maxCutoff * maxCutoff;
        auto d2 = ctx->getDistSquaredFun();

        for (auto &cell : cells) {
            cell.particleIndices.clear();
        }
        neighbor_list.clear();
        neighbor_list.resize(data.size());
        data_t::index_t idx = 0;

        for (const auto &entry : data) {
            auto cell = getCell(entry.pos);
            if (cell) {
                cell->particleIndices.push_back(idx);
            }
            ++idx;
        }

        using cell_it_t = decltype(cells.begin());
        const auto size = cells.size();
        const std::size_t grainSize = size / config->nThreads();

        auto worker = [this, c2, d2](cell_it_t cellsBegin, cell_it_t cellsEnd) -> void {
            for (cell_it_t _b = cellsBegin; _b != cellsEnd; ++_b) {
                auto &cell = *_b;
                setUpCell(cell, c2, d2);
            }
        };

        std::vector<thd::scoped_thread> threads;
        threads.reserve(config->nThreads());

        auto it_cells = cells.begin();
        for (int i = 0; i < config->nThreads() - 1; ++i) {
            threads.push_back(thd::scoped_thread(std::thread(worker, it_cells, it_cells + grainSize)));
            it_cells += grainSize;
        }
        threads.push_back(thd::scoped_thread(std::thread(worker, it_cells, cells.end())));
    }
}

void CPUDNeighborList::create() {
    simBoxSize = ctx->getBoxSize();
    setupCells();
    fillCells();
}

CPUDNeighborList::Cell *CPUDNeighborList::getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k) {
    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
    else if (i < 0 || i >= nCells[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
    else if (j < 0 || j >= nCells[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
    else if (k < 0 || k >= nCells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, nCells[1], nCells[2]);
    if (cix < cells.size()) {
        return &cells.at(static_cast<cell_index>(cix));
    } else {
        log::critical("CPUDNeighborList::getCell(nonconst): Requested cell ({},{},{})={}, but there are "
                                         "only {} cells.", i, j, k, cix, cells.size());
        throw std::runtime_error("tried to get cell index that was too large");
    }
}

CPUDNeighborList::Cell *CPUDNeighborList::getCell(const readdy::model::Particle::pos_type &pos) {
    const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    return getCell(i, j, k);
}

const std::vector<CPUDNeighborList::neighbor_t> &CPUDNeighborList::neighbors(CPUDNeighborList::particle_index const entry) const {
    if (maxCutoff > 0) {
        return neighbor_list.at(entry);
    }
    return no_neighbors;
}

const CPUDNeighborList::Cell *const CPUDNeighborList::getCell(const readdy::model::Particle::pos_type &pos) const {
    const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    return getCell(i, j, k);
}

const CPUDNeighborList::Cell *const CPUDNeighborList::getCell(CPUDNeighborList::signed_cell_index i,
                                                      CPUDNeighborList::signed_cell_index j,
                                                      CPUDNeighborList::signed_cell_index k) const {

    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
    else if (i < 0 || i >= nCells[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
    else if (j < 0 || j >= nCells[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
    else if (k < 0 || k >= nCells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, nCells[1], nCells[2]);
    if (cix < cells.size()) {
        return &cells.at(static_cast<cell_index>(cix));
    } else {
        log::critical("CPUDNeighborList::getCell(const): Requested cell ({},{},{})={}, but there are "
                                         "only {} cells.", i, j, k, cix, cells.size());
        throw std::out_of_range("tried to access an invalid cell");
    }
}

const std::vector<CPUDNeighborList::neighbor_t> &CPUDNeighborList::find_neighbors(particle_index const entry) const {
    if (maxCutoff > 0 && entry < neighbor_list.size()) {
        return neighbor_list.at(entry);
    }
    return no_neighbors;
}

CPUDNeighborList::~CPUDNeighborList() = default;

void CPUDNeighborList::Cell::addNeighbor(CPUDNeighborList::Cell *cell) {
    if (cell && cell->contiguous_index != contiguous_index
        && (enoughCells || std::find(neighbors.begin(), neighbors.end(), cell) == neighbors.end())) {
        neighbors.push_back(cell);
    }
}

CPUDNeighborList::Cell::Cell(cell_index i, cell_index j, cell_index k,
                         const std::array<CPUDNeighborList::cell_index, 3> &nCells)
        : contiguous_index(get_contiguous_index(i, j, k, nCells[1], nCells[2])),
          enoughCells(nCells[0] >= 5 && nCells[1] >= 5 && nCells[2] >= 5) {
}

bool operator==(const CPUDNeighborList::Cell &lhs, const CPUDNeighborList::Cell &rhs) {
    return lhs.contiguous_index == rhs.contiguous_index;
}

bool operator!=(const CPUDNeighborList::Cell &lhs, const CPUDNeighborList::Cell &rhs) {
    return !(lhs == rhs);
}

void CPUDNeighborList::setUpCell(CPUDNeighborList::Cell &cell, const double cutoffSquared, const ctx_t::dist_squared_fun &d2) {
    for (const auto &pI : cell.particleIndices) {
        auto &entry_i = data.entry_at(pI);
        auto &neighbors_i = neighbor_list.at(pI);
        neighbors_i.clear();
        for (const auto &pJ : cell.particleIndices) {
            if (pI != pJ) {
                const auto distSquared = d2(entry_i.pos, data.centry_at(pJ).pos);
                if (distSquared < cutoffSquared) {
                    neighbors_i.push_back({pJ, distSquared});
                }
            }
        }
        for (const auto &neighboringCell : cell.neighbors) {
            for (const auto pJ : neighboringCell->particleIndices) {
                const auto distSquared = d2(entry_i.pos, data.centry_at(pJ).pos);
                if (distSquared < cutoffSquared) {
                    neighbors_i.push_back({pJ, distSquared});
                }
            }
        }

    }
}

int CPUDNeighborList::getCellIndex(const readdy::model::Particle::pos_type &pos) const {
    signed_cell_index i = static_cast<const signed_cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    signed_cell_index j = static_cast<const signed_cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    signed_cell_index k = static_cast<const signed_cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
    else if (i < 0 || i >= nCells[0]) return -1;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
    else if (j < 0 || j >= nCells[1]) return -1;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
    else if (k < 0 || k >= nCells[2]) return -1;
    return get_contiguous_index(i, j, k, nCells[1], nCells[2]);
}

const double CPUDNeighborList::getMaxCutoff() const {
    return maxCutoff;
}

CPUDNeighborList::iterator CPUDNeighborList::begin() {
        return neighbor_list.begin();
}

CPUDNeighborList::iterator CPUDNeighborList::end() {
    return neighbor_list.end();
}

CPUDNeighborList::const_iterator CPUDNeighborList::cbegin() const {
    return neighbor_list.cbegin();
}

CPUDNeighborList::const_iterator CPUDNeighborList::cend() const {
    return neighbor_list.cend();
}

CPUDNeighborList::const_iterator CPUDNeighborList::begin() const {
    return neighbor_list.cbegin();
}

CPUDNeighborList::const_iterator CPUDNeighborList::end() const {
    return neighbor_list.cend();
}

}
}
}
}