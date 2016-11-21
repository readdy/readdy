/**
 * << detailed description >>
 *
 * @file NeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.07.16
 */

#ifndef READDY_CPUKERNEL_NEIGHBORLIST_H
#define READDY_CPUKERNEL_NEIGHBORLIST_H

#include <memory>

#include <readdy/common/make_unique.h>
#include <readdy/model/KernelContext.h>

#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/kernel/cpu/util/Config.h>

#include "ParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

class NeighborList {
public: 
    class Cell;
    using cell_index = unsigned int;
    using signed_cell_index = typename std::make_signed<cell_index>::type;
    using particle_index = ParticleData::index_t;
    using neighbor_t = ParticleData::Neighbor;
    using data_t = readdy::kernel::cpu::model::ParticleData;
    using ctx_t = readdy::model::KernelContext;
    using container_t = std::vector<std::vector<neighbor_t>>;
    using skin_size_t = double;
    using data_iter_t = decltype(std::declval<data_t>().entries.begin());
    using hilbert_index_t = unsigned int;
    using cell_iter_t = decltype(std::declval<std::vector<Cell>>().begin());
private: 
    std::vector<Cell> cells;
    skin_size_t skin_size;
    bool initialSetup = true;
public:
    struct Cell {
        std::vector<Cell *> neighbors{};
        std::vector<particle_index> particleIndices{};
        double maximal_displacements[2];
        cell_index contiguous_index;
        bool enoughCells;

        // dirty flag indicating whether the cell and its neighboring cells have to be re-created
        bool dirty {false};

        Cell(cell_index i, cell_index j, cell_index k, const std::array<cell_index, 3> &nCells);

        void addNeighbor(Cell *cell);

        void checkDirty(skin_size_t skin);

        friend bool operator==(const Cell &lhs, const Cell &rhs);

        friend bool operator!=(const Cell &lhs, const Cell &rhs);
    };

    NeighborList(const ctx_t *const context, data_t &data, util::Config const *const config, skin_size_t = 0);

    virtual ~NeighborList();

    virtual void setupCells();

    virtual void setupNeighboringCells(const signed_cell_index i, const signed_cell_index j, const signed_cell_index k);

    void clear();

    const std::vector<neighbor_t>& neighbors(const particle_index entry) const;
    const std::vector<neighbor_t>& find_neighbors(const particle_index) const;

    virtual void fillCells();

    virtual void create();

    void updateData(data_t::update_t&& update);

    void displace(data_iter_t iter, const readdy::model::Vec3& vec);
    void displace(data_t::Entry&, const data_t::particle_type::pos_type& delta);
    void displace(data_t::index_t entry, const data_t::particle_type::pos_type& delta);
    void setPosition(data_t::index_t entry, data_t::particle_type::pos_type&& newPosition);

    void remove(const particle_index);

    void insert(const particle_index);

    auto begin() -> decltype(std::declval<data_t>().neighbors.begin()) {
        return data.neighbors.begin();
    }
    auto end() -> decltype(std::declval<data_t>().neighbors.end()) {
        return data.neighbors.end();
    }
    auto cbegin() const -> decltype(std::declval<data_t>().neighbors.cbegin()) {
        return data.neighbors.cbegin();
    }
    auto cend() const -> decltype(std::declval<data_t>().neighbors.cend()) {
        return data.neighbors.cend();
    }
    auto begin() const -> decltype(std::declval<data_t>().neighbors.cbegin()) {
        return cbegin();
    }
    auto end() const -> decltype(std::declval<data_t>().neighbors.cend()) {
        return cend();
    }

    void setSkinSize(skin_size_t skin_size);

    hilbert_index_t getHilbertIndex(std::size_t i, std::size_t j, std::size_t k) const;

protected:

    const readdy::model::KernelContext *const ctx;

    using cell_size_t = decltype(ctx->getBoxSize());

    cell_size_t simBoxSize;

    std::array<cell_index, 3> nCells{{0, 0, 0}};
    readdy::model::Vec3 cellSize{0, 0, 0};
    double maxCutoff = 0;
    double maxCutoffPlusSkin = 0;
    util::Config const *const config;

    const Cell* const getCell(const readdy::model::Particle::pos_type &pos) const;

    int getCellIndex(const readdy::model::Particle::pos_type &pos) const;

    Cell *getCell(const readdy::model::Particle::pos_type &pos);
    Cell *getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k);

    const Cell * const getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k) const;

    std::unordered_set<Cell*> findDirtyCells();

    void setUpCell(NeighborList::Cell &cell, const double cutoffSquared, const ctx_t::dist_squared_fun& d2);

    data_t& data;
};


}
}
}
}
#endif //READDY_CPUKERNEL_NEIGHBORLIST_H
