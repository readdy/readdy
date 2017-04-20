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
 * @date 14.07.16
 */

#pragma once
#include <memory>

#include <readdy/common/make_unique.h>
#include <readdy/model/KernelContext.h>

#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/common/thread/Config.h>
#include <readdy/common/signals.h>

#include "CPUParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

class CPUNeighborList {
public:

    using reorder_signal_t = readdy::signals::signal<void(const std::vector<std::size_t>)>;

    struct Cell;
    struct CoarseCell;
    using cell_index = unsigned int;
    using signed_cell_index = typename std::make_signed<cell_index>::type;
    using particle_index = CPUParticleData::index_t;
    using neighbor_t = CPUParticleData::Neighbor;
    using data_t = readdy::kernel::cpu::model::CPUParticleData;
    using ctx_t = readdy::model::KernelContext;
    using container_t = std::vector<std::vector<neighbor_t>>;
    using skin_size_t = double;
    using data_iter_t = data_t::iterator;
    using hilbert_index_t = unsigned int;
    using cell_iter_t = decltype(std::declval<std::vector<Cell>>().begin());

    using iterator = decltype(std::declval<data_t>().neighbors.begin());
    using const_iterator = decltype(std::declval<data_t>().neighbors.cbegin());

private:
    std::vector<Cell> cells;
    std::vector<CoarseCell> coarseCells;
    skin_size_t skin_size;
    bool initialSetup = true;
public:

    CPUNeighborList(const ctx_t *const context, data_t &data, readdy::util::thread::Config const *const config,
                    skin_size_t = 0);

    virtual ~CPUNeighborList();

    virtual void setupCells();

    virtual void setupNeighboringCells(const signed_cell_index i, const signed_cell_index j, const signed_cell_index k);

    void clear();

    const std::vector<neighbor_t> &neighbors(const particle_index entry) const;

    const std::vector<neighbor_t> &find_neighbors(const particle_index) const;

    virtual void fillCells();

    virtual void create();

    void updateData(data_t::update_t &&update);

    void displace(data_iter_t iter, const readdy::model::Vec3 &vec);

    void displace(data_t::Entry &, const data_t::particle_type::pos_type &delta);

    void displace(data_t::index_t entry, const data_t::particle_type::pos_type &delta);

    void setPosition(data_t::index_t entry, data_t::particle_type::pos_type &&newPosition);

    readdy::signals::scoped_connection registerReorderEventListener(const reorder_signal_t::slot_type &slot);

    void remove(const particle_index);

    void insert(const particle_index);

    iterator begin() {
        return data.neighbors.begin();
    }

    iterator end() {
        return data.neighbors.end();
    }

    const_iterator cbegin() const {
        return data.neighbors.cbegin();
    }

    const_iterator cend() const {
        return data.neighbors.cend();
    }

    const_iterator begin() const {
        return cbegin();
    }

    const_iterator end() const {
        return cend();
    }

    void setSkinSize(skin_size_t skin_size);

    void setGroupParticlesOnCreation(bool groupParticlesOnCreation);

    hilbert_index_t getHilbertIndex(std::size_t i, std::size_t j, std::size_t k) const;

    std::unordered_set<Cell *> findDirtyCells();

    void handleDirtyCells(std::unordered_set<Cell*> dirties);

    std::size_t hash_pos(const data_t::particle_type::pos_type &pos) const;

    std::tuple<cell_index, cell_index, cell_index> mapPositionToCell(const data_t::particle_type::pos_type &pos) const;

protected:

    bool isInCell(const Cell &cell, const data_t::particle_type::pos_type &pos) const;

    const readdy::model::KernelContext *const ctx;

    using cell_size_t = decltype(ctx->getBoxSize());

    cell_size_t simBoxSize;

    std::array<cell_index, 3> nCells{{0, 0, 0}};
    readdy::model::Vec3 cellSize{0, 0, 0};
    double maxCutoff = 0;
    double maxCutoffPlusSkin = 0;
    readdy::util::thread::Config const *const config;

    const Cell *const getCell(const readdy::model::Particle::pos_type &pos) const;

    int getCellIndex(const readdy::model::Particle::pos_type &pos) const;

    Cell *getCell(const readdy::model::Particle::pos_type &pos);

    Cell *getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k);

    const Cell *const getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k) const;

    void setUpCell(CPUNeighborList::Cell &cell, const double cutoffSquared, const ctx_t::dist_squared_fun &d2);

    data_t &data;

    std::unique_ptr<reorder_signal_t> reorderSignal;

    bool groupParticlesOnCreation;
};


}
}
}
}
