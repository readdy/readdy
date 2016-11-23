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

#include "ParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace model {

class NeighborList {
public: 
    class Cell;
    class Neighbor;
    using cell_index = unsigned int;
    using signed_cell_index = typename std::make_signed<cell_index>::type;
    using particle_index = ParticleData::index_t;
    using neighbor_t = Neighbor;
    using data_t = readdy::kernel::cpu_dense::model::ParticleData;
    using ctx_t = readdy::model::KernelContext;
    using container_t = std::vector<std::vector<neighbor_t>>;
    using data_iter_t = decltype(std::declval<data_t>().begin());
    using cell_iter_t = decltype(std::declval<std::vector<Cell>>().begin());

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
    std::vector<std::vector<Neighbor>> neighbor_list;
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

    NeighborList(const ctx_t *const, data_t &, readdy::util::thread::Config const *const);

    virtual ~NeighborList();

    virtual void setupCells();

    virtual void setupNeighboringCells(const signed_cell_index i, const signed_cell_index j, const signed_cell_index k);

    void clear();

    const std::vector<neighbor_t>& neighbors(const particle_index entry) const;
    const std::vector<neighbor_t>& find_neighbors(const particle_index) const;

    virtual void fillCells();

    virtual void create();

    auto begin() -> decltype(neighbor_list.begin()) {
        return neighbor_list.begin();
    }
    auto end() -> decltype(neighbor_list.end()) {
        return neighbor_list.end();
    }
    auto cbegin() const -> decltype(neighbor_list.cbegin()) {
        return neighbor_list.cbegin();
    }
    auto cend() const -> decltype(neighbor_list.cend()) {
        return neighbor_list.cend();
    }
    auto begin() const -> decltype(neighbor_list.cbegin()) {
        return cbegin();
    }
    auto end() const -> decltype(neighbor_list.cend()) {
        return cend();
    }

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

    void setUpCell(NeighborList::Cell &cell, const double cutoffSquared, const ctx_t::dist_squared_fun& d2);

    data_t& data;
};


}
}
}
}
#endif //READDY_DENSE_NEIGHBORLIST_H
