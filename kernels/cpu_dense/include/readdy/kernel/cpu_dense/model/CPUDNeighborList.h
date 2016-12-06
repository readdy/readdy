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
 * @file NeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_DENSE_NEIGHBORLIST_H
#define READDY_DENSE_NEIGHBORLIST_H

#include <memory>

#include <readdy/common/make_unique.h>
#include <readdy/model/KernelContext.h>

#include <readdy/kernel/cpu_dense/model/ParticleIndexPair.h>
#include <readdy/common/thread/Config.h>

#include "CPUDParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace model {

class CPUDNeighborList {
public:
    class Cell;
    class Neighbor;
    using cell_index = unsigned int;
    using signed_cell_index = typename std::make_signed<cell_index>::type;
    using particle_index = CPUDParticleData::index_t;
    using neighbor_t = Neighbor;
    using data_t = readdy::kernel::cpu_dense::model::CPUDParticleData;
    using ctx_t = readdy::model::KernelContext;
    using container_t = std::vector<std::vector<neighbor_t>>;
    using data_iter_t = decltype(std::declval<data_t>().begin());
    using cell_iter_t = decltype(std::declval<std::vector<Cell>>().begin());
    using neighbor_list_t = std::vector<std::vector<Neighbor>>;

    using iterator = neighbor_list_t::iterator;
    using const_iterator = neighbor_list_t::const_iterator;

    struct Neighbor {
        using index_t = data_t::entries_t::size_type;
        index_t idx;
        double d2;

        Neighbor(const index_t idx, const double d2);
        Neighbor(const Neighbor&) = delete;
        Neighbor& operator=(const Neighbor&) = delete;
        Neighbor(Neighbor&&) = default;
        Neighbor& operator=(Neighbor&&) = delete;
    };


private: 
    std::vector<Cell> cells;
    neighbor_list_t neighbor_list;
public:
    struct Cell {
        std::vector<Cell *> neighbors{};
        std::vector<particle_index> particleIndices{};
        cell_index contiguous_index;
        bool enoughCells;

        Cell(cell_index i, cell_index j, cell_index k, const std::array<cell_index, 3> &nCells);

        void addNeighbor(Cell *cell);

        friend bool operator==(const Cell &lhs, const Cell &rhs);

        friend bool operator!=(const Cell &lhs, const Cell &rhs);
    };

    CPUDNeighborList(const ctx_t *const, data_t &, readdy::util::thread::Config const *const);

    virtual ~CPUDNeighborList();

    virtual void setupCells();

    virtual void setupNeighboringCells(const signed_cell_index i, const signed_cell_index j, const signed_cell_index k);

    void clear();

    const std::vector<neighbor_t>& neighbors(const particle_index entry) const;
    const std::vector<neighbor_t>& find_neighbors(const particle_index) const;

    virtual void fillCells();

    virtual void create();

    iterator begin();
    iterator end();
    const_iterator cbegin() const;
    const_iterator cend() const;
    const_iterator begin() const;
    const_iterator end() const;

    const double getMaxCutoff() const;

protected:

    const readdy::model::KernelContext *const ctx;

    using cell_size_t = decltype(ctx->getBoxSize());

    cell_size_t simBoxSize;

    std::array<cell_index, 3> nCells{{0, 0, 0}};
    readdy::model::Vec3 cellSize{0, 0, 0};
    double maxCutoff = 0;
    readdy::util::thread::Config const *const config;

    const Cell* const getCell(const readdy::model::Particle::pos_type &pos) const;

    int getCellIndex(const readdy::model::Particle::pos_type &pos) const;

    Cell *getCell(const readdy::model::Particle::pos_type &pos);
    Cell *getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k);

    const Cell * const getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k) const;

    void setUpCell(CPUDNeighborList::Cell &cell, const double cutoffSquared, const ctx_t::dist_squared_fun& d2);

    data_t& data;
};


}
}
}
}
#endif //READDY_DENSE_NEIGHBORLIST_H
