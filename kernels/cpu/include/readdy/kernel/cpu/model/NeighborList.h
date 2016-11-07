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
#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>
#include <readdy/kernel/cpu/util/scoped_thread.h>
#include <readdy/kernel/cpu/util/Config.h>
#include "ParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {
struct Neighbor {
    using index_t = readdy::kernel::cpu::model::ParticleData::Entry *;
    index_t idx;
    double d2;

    Neighbor(const index_t idx, const double d2);
    Neighbor(const Neighbor&) = delete;
    Neighbor& operator=(const Neighbor&) = delete;
    Neighbor(Neighbor&&);
    Neighbor& operator=(Neighbor&&);
};

class NeighborList {
public: 
    class Cell;
    using cell_index = unsigned int;
    using signed_cell_index = typename std::make_signed<cell_index>::type;
    using particle_index = Neighbor::index_t;
    using neighbor_t = Neighbor;
    using data_t = readdy::kernel::cpu::model::ParticleData;
    using container_t = std::unordered_map<particle_index, std::vector<neighbor_t>>;
private: 
    std::vector<Cell> cells;
    std::vector<container_t> maps;
public:
    struct Cell {
        std::vector<Cell *> neighbors{};
        std::vector<particle_index> particleIndices{};
        const cell_index id;
        const bool enoughCells;

        // dirty flag indicating whether the cell and its neighboring cells have to be re-created
        bool dirty = true;

        Cell(cell_index i, cell_index j, cell_index k, const std::array<cell_index, 3> &nCells);

        void addNeighbor(Cell *cell);

        friend bool operator==(const Cell &lhs, const Cell &rhs) {
            return lhs.id == rhs.id;
        }

        friend bool operator!=(const Cell &lhs, const Cell &rhs) {
            return !(lhs == rhs);
        }
    };

    NeighborList(const readdy::model::KernelContext *const context, util::Config const *const config);

    virtual ~NeighborList();

    virtual void setupCells();

    virtual void setupNeighboringCells(const signed_cell_index i, const signed_cell_index j, const signed_cell_index k);

    void clear();

    const std::vector<neighbor_t>& neighbors(const particle_index entry) const;
    const std::vector<neighbor_t>& find_neighbors(const particle_index) const;

    virtual void fillCells(data_t &data);

    virtual void create(data_t &data);

    void updateData(data_t &data, data_t::update_t update);

    void remove(const particle_index);

    void insert(const data_t &data, const particle_index);

    auto begin() -> decltype(maps.begin()) {
        return maps.begin();
    }
    auto end() -> decltype(maps.end()) {
        return maps.end();
    }
    auto cbegin() const -> decltype(maps.cbegin()) {
        return maps.cbegin();
    }
    auto cend() const -> decltype(maps.cend()) {
        return maps.cend();
    }
    auto begin() const -> decltype(maps.cbegin()) {
        return cbegin();
    }
    auto end() const -> decltype(maps.cend()) {
        return cend();
    }
    auto n_maps() const -> decltype(maps.size()) {
        return maps.size();
    }

protected:

    const readdy::model::KernelContext *const ctx;

    using cell_size_t = decltype(ctx->getBoxSize());

    cell_size_t simBoxSize;

    std::array<cell_index, 3> nCells{{0, 0, 0}};
    readdy::model::Vec3 cellSize{0, 0, 0};
    double maxCutoff = 0;
    util::Config const *const config;

    const Cell* const getCell(const readdy::model::Particle::pos_type &pos) const;

    Cell *getCell(const readdy::model::Particle::pos_type &pos);

    Cell *getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k);

    const Cell * const getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k) const;

    util::Config::n_threads_t getMapsIndex(const Cell* const cell) const;

    container_t& getPairs(const Cell* const cell);

    const container_t& getPairs(const Cell* const cell) const;
};
}
}
}
}
#endif //READDY_CPUKERNEL_NEIGHBORLIST_H
