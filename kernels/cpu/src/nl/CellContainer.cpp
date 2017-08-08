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
 * @file CellContainer.cpp
 * @brief 
 * @author clonker
 * @date 4/21/17
 */

#include <cmath>
#include <atomic>

#include <readdy/common/numeric.h>
#include <readdy/kernel/cpu/nl/CellContainer.h>
#include <readdy/kernel/cpu/nl/SubCell.h>
#include <readdy/kernel/cpu/util/config.h>
#include <readdy/common/thread/barrier.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

CellContainer::sub_cells_t &CellContainer::sub_cells() {
    return _sub_cells;
}

const CellContainer::sub_cells_t &CellContainer::sub_cells() const {
    return _sub_cells;
}

CellContainer::dimension &CellContainer::n_sub_cells() {
    return _n_sub_cells;
}

const CellContainer::dimension &CellContainer::n_sub_cells() const {
    return _n_sub_cells;
}

CellContainer::vec3 &CellContainer::size() {
    return _size;
}

const CellContainer::vec3 &CellContainer::size() const {
    return _size;
}

CellContainer::cell_ref_list &CellContainer::neighbors() {
    return _neighbors;
}

const CellContainer::cell_ref_list &CellContainer::neighbors() const {
    return _neighbors;
}

CellContainer::displacement_arr &CellContainer::maximal_displacements() {
    return _maximal_displacements;
}

const CellContainer::displacement_arr &CellContainer::maximal_displacements() const {
    return _maximal_displacements;
}

CellContainer::cell_index &CellContainer::contiguous_index() {
    return _contiguous_index;
}

const CellContainer::cell_index &CellContainer::contiguous_index() const {
    return _contiguous_index;
}

CellContainer::CellContainer(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                             const readdy::util::thread::Config &config)
        : _data(data), _context(context), _config(config), _size(context.getBoxSize()),
          _root_size(context.getBoxSize()) {}

CellContainer::sub_cell *const CellContainer::sub_cell_for_position(const CellContainer::vec3 &pos) {
    return const_cast<CellContainer::sub_cell *const>(static_cast<const CellContainer *>(this)->sub_cell_for_position(
            pos));
}

const CellContainer::sub_cell *const CellContainer::sub_cell_for_position(const CellContainer::vec3 &pos) const {
    const auto i = static_cast<const cell_index>(floor(
            (pos.x + .5 * _root_size.x - _offset.x) / _sub_cell_size.x));
    const auto j = static_cast<const cell_index>(floor(
            (pos.y + .5 * _root_size.y - _offset.y) / _sub_cell_size.y));
    const auto k = static_cast<const cell_index>(floor(
            (pos.z + .5 * _root_size.z - _offset.z) / _sub_cell_size.z));
    return sub_cell_for_grid_index(std::tie(i, j, k));
}

CellContainer::sub_cell *const CellContainer::sub_cell_for_grid_index(const CellContainer::grid_index &index) {
    return const_cast<CellContainer::sub_cell *const>(
            static_cast<const CellContainer *>(this)->sub_cell_for_grid_index(index)
    );
}

const CellContainer::sub_cell *const
CellContainer::sub_cell_for_grid_index(const CellContainer::grid_index &index) const {
    auto i = std::get<0>(index);
    auto j = std::get<1>(index);
    auto k = std::get<2>(index);
    if (i < 0 || i >= _n_sub_cells[0]) return nullptr;
    if (j < 0 || j >= _n_sub_cells[1]) return nullptr;
    if (k < 0 || k >= _n_sub_cells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, _n_sub_cells[1], _n_sub_cells[2]);
    if (cix < _sub_cells.size()) {
        return &_sub_cells.at(static_cast<cell_index>(cix));
    }
    log::critical("CPUNeighborList::getCell(const): Requested cell ({},{},{})={}, but there are "
                              "only {} cells.", i, j, k, cix, _sub_cells.size());
    throw std::runtime_error("tried to get cell index that was too large");
}

void CellContainer::update_displacements() {
    update_sub_cell_displacements();
    _maximal_displacements[0] = 0;
    _maximal_displacements[1] = 0;
    for (auto &cell : _sub_cells) {
        const auto &updated_displacements = cell.maximal_displacements();
        if (updated_displacements[0] > _maximal_displacements[0]) {
            _maximal_displacements[0] = updated_displacements[0];
            _maximal_displacements[1] = std::max(_maximal_displacements[1], updated_displacements[1]);
        } else if (updated_displacements[0] > _maximal_displacements[1]) {
            _maximal_displacements[1] = updated_displacements[0];
        }
    }
}

void CellContainer::subdivide(const scalar desired_cell_width) {
    for (std::uint8_t i = 0; i < 3; ++i) {
        _n_sub_cells[i] = std::max(static_cast<cell_index>(1),
                                   static_cast<cell_index>(std::floor(_size[i] / desired_cell_width)));
        _sub_cell_size[i] = _size[i] / static_cast<scalar>(_n_sub_cells[i]);
    }
    log::debug("resulting cell size = {}", _sub_cell_size);
    _sub_cells.reserve(_n_sub_cells[0] * _n_sub_cells[1] * _n_sub_cells[2]);
    for (cell_index i = 0; i < _n_sub_cells[0]; ++i) {
        for (cell_index j = 0; j < _n_sub_cells[1]; ++j) {
            for (cell_index k = 0; k < _n_sub_cells[2]; ++k) {
                _sub_cells.emplace_back(this, vec3{i * _sub_cell_size.x, j * _sub_cell_size.y, k * _sub_cell_size.z});
                _sub_cells.back()._size = _sub_cell_size;
                _sub_cells.back()._contiguous_index = get_contiguous_index(i, j, k, _n_sub_cells[1], _n_sub_cells[2]);
            }
        }
    }
}

const model::CPUParticleData &CellContainer::data() const {
    return _data;
}

const readdy::model::KernelContext &CellContainer::context() const {
    return _context;
}

void CellContainer::update_sub_cell_displacements() {
    execute_for_each_sub_cell([](sub_cell &cell) {
        cell.update_displacements();
    });
}

const readdy::util::thread::Config &CellContainer::config() const {
    return _config;
}

const CellContainer::level_t CellContainer::level() const {
    return _level;
}

const CellContainer *const CellContainer::root() const {
    auto root = this;
    while (root->_super_cell != nullptr) root = root->_super_cell;
    return root;
}

CellContainer *const CellContainer::root() {
    auto root = this;
    while (root->_super_cell != nullptr) root = root->_super_cell;
    return root;
}

void CellContainer::setup_uniform_neighbors() {
    execute_for_each_sub_cell([](sub_cell &cell) {
        cell.setup_uniform_neighbors(1);
    });
}

void CellContainer::refine_uniformly() {
    execute_for_each_sub_cell([](sub_cell &cell) {
        cell.refine_uniformly();
    });
}

CellContainer::sub_cell *const CellContainer::leaf_cell_for_position(const CellContainer::vec3 &pos) {
    return const_cast<CellContainer::sub_cell *const>(
            static_cast<const CellContainer *>(this)->leaf_cell_for_position(pos)
    );
}

const CellContainer::sub_cell *const CellContainer::leaf_cell_for_position(const CellContainer::vec3 &pos) const {
    auto cell = sub_cell_for_position(pos);
    if (cell != nullptr) {
        while (!cell->is_leaf()) cell = cell->sub_cell_for_position(pos);
    }
    return cell;
}

CellContainer::sub_cell *const
CellContainer::sub_cell_for_position(const vec3 &pos, const level_t level) {
    return const_cast<CellContainer::sub_cell *const>(static_cast<const CellContainer *>(this)->sub_cell_for_position(
            pos, level));
}

const CellContainer::sub_cell *const
CellContainer::sub_cell_for_position(const vec3 &pos, const level_t level) const {
    if (level == 0) throw std::invalid_argument("level needs to be larger 0");
    auto cell = sub_cell_for_position(pos);
    if (cell != nullptr) {
        while (cell->level() != level) cell = cell->sub_cell_for_position(pos);
    }
    return cell;
}

const CellContainer::vec3 &CellContainer::offset() const {
    return _offset;
}

const CellContainer *const CellContainer::super_cell() const {
    return _super_cell;
}

void CellContainer::insert_particle(const CellContainer::particle_index index, bool mark_dirty) const {
    const auto &entry = data().entry_at(index);
    if (!entry.is_deactivated()) {
        auto cell = leaf_cell_for_position(entry.position());
        if (cell != nullptr) {
            cell->insert_particle(index);
            if(mark_dirty) {
                cell->super_cell()->set_dirty();
            }
        }
    } else {
        log::critical("Tried inserting a deactivated particle ({}) into the neighbor list!", index);
    }
}

void CellContainer::insert_particle(const CellContainer::particle_index index, bool mark_dirty) {
    const auto &entry = data().entry_at(index);
    if (!entry.is_deactivated()) {
        auto cell = leaf_cell_for_position(entry.position());
        if (cell != nullptr) {
            cell->insert_particle(index);
            if(mark_dirty) {
                cell->super_cell()->set_dirty();
            }
        }
    } else {
        log::critical("Tried inserting a deactivated particle ({}) into the neighbor list!", index);
    }
}

void CellContainer::clear() {
    reset_max_displacements();
    execute_for_each_sub_cell([](sub_cell &cell) {
        cell.clear();
    });
}

bool CellContainer::update_sub_cell_displacements_and_mark_dirty(const scalar cutoff, const scalar skin) {
    update_sub_cell_displacements();
    std::atomic<bool> status{true};
    std::atomic<std::size_t> n_dirty{0};
    execute_for_each_sub_cell([&](sub_cell &cell) {
        displacement_arr max_max_displacements (cell.maximal_displacements());
        for(const auto& neighbor_cell : cell.neighbors()) {
            if(max_max_displacements[0] < neighbor_cell->maximal_displacements()[0]) {
                max_max_displacements[0] = neighbor_cell->maximal_displacements()[0];
                if(max_max_displacements[1] < neighbor_cell->maximal_displacements()[1]) {
                    max_max_displacements[1] = neighbor_cell->maximal_displacements()[1];
                }
            } else if(max_max_displacements[1] < neighbor_cell->maximal_displacements()[0]) {
                max_max_displacements[1] = neighbor_cell->maximal_displacements()[0];
            }
        }
        if ((max_max_displacements[0] + max_max_displacements[1]) >= skin) {
            if (max_max_displacements[0] > skin + cutoff) {
                status = false;
            }
            cell.set_dirty();
            //cell.reset_particles_displacements();
            /*for (const auto &neighbor_cell : cell.neighbors()) {
                neighbor_cell->set_dirty();
            }*/
        }
    });
    execute_for_each_sub_cell([&](sub_cell &cell) {
        if (cell.is_dirty()) {
            std::atomic_fetch_add<std::size_t>(&n_dirty, 1);
            cell.reset_max_displacements();
        }
    });
    _n_dirty_macro_cells = n_dirty.load();
    return status.load();
}

void CellContainer::update_dirty_cells() {
    using barrier_t = util::thread::barrier;
    // todo go through macro cells, see if dirty => the current particles might be in the neighboring cells

    // 1. reset displacement of DIRTY cells
    // 3. reassign these particles to cells, rebuild verlet list (in NL?) for particles in DIRTY and NEIGHBORS

    auto worker = [this](std::size_t, const CellContainer& container,
                         const barrier_t& barrier,
                         sub_cells_t::iterator begin, sub_cells_t::iterator end) {
        std::vector<std::vector<particle_index>> dirty_cell_particles;
        {
            for (auto it = begin; it != end; ++it) {
                if (it->is_dirty()) {
                    it->reset_particles_displacements();
                    dirty_cell_particles.emplace_back(it->collect_contained_particles());
                    it->clear();
                }
            }
        }
        barrier.wait();
        {
            for (const auto &dirty_cell : dirty_cell_particles) {
                for (const auto p_idx : dirty_cell) {
                    auto &entry = _data.entry_at(p_idx);
                    // entry.displacement = 0;
                    auto cell = container.leaf_cell_for_position(entry.position());
                    if (cell != nullptr) {
                        cell->insert_particle(p_idx);
                    }
                }
            }
        }
    };

    barrier_t barrier{_config.nThreads()};
    const auto grainSize = sub_cells().size() / config().nThreads();

    std::vector<std::function<void(std::size_t)>> executables;
    executables.reserve(config().nThreads());
    const auto& executor = *config().executor();

    auto it = sub_cells().begin();
    for (auto i = 0; i < config().nThreads() - 1; ++i) {
        executables.push_back(executor.pack(worker, std::cref(*this), std::cref(barrier), it, it + grainSize));
        it += grainSize;
    }
    executables.push_back(executor.pack(worker, std::cref(*this), std::cref(barrier), it, sub_cells().end()));
    executor.execute_and_wait(std::move(executables));
}

const std::size_t CellContainer::n_dirty_macro_cells() const {
    return _n_dirty_macro_cells;
}

const std::size_t CellContainer::n_sub_cells_total() const {
    return _n_sub_cells[0] * _n_sub_cells[1] * _n_sub_cells[2];
}

void CellContainer::execute_for_each_leaf(const std::function<void(const CellContainer::sub_cell &)> &function) {
    execute_for_each_sub_cell([&](sub_cell &cell) {
        for (const auto &sub_cell : cell.sub_cells()) {
            function(sub_cell);
        }
    });
}

void CellContainer::execute_for_each_sub_cell(const std::function<void(CellContainer::sub_cell &)> &function) {
    const auto grainSize = sub_cells().size() / config().nThreads();
    auto worker = [&](std::size_t, CellContainer::sub_cells_t::iterator begin, CellContainer::sub_cells_t::iterator end) {
        for (auto it = begin; it != end; ++it) {
            function(*it);
        }
    };

    const auto& executor = *config().executor();
    std::vector<std::function<void(std::size_t)>> executables;
    executables.reserve(config().nThreads());

    auto it = sub_cells().begin();
    for (auto i = 0; i < config().nThreads() - 1; ++i) {
        executables.push_back(executor.pack(worker, it, it+grainSize));
        it += grainSize;
    }
    executables.push_back(executor.pack(worker, it, sub_cells().end()));
    executor.execute_and_wait(std::move(executables));
}

model::CPUParticleData &CellContainer::data() {
    return _data;
}

void CellContainer::unset_dirty() const {
    execute_for_each_sub_cell([](const sub_cell &cell) {
        cell.unset_dirty();
    });
}

void CellContainer::reset_max_displacements() {
    maximal_displacements()[0] = 0;
    maximal_displacements()[1] = 0;
    execute_for_each_sub_cell([](sub_cell &cell) {
        cell.reset_max_displacements();
    });
}

void CellContainer::set_dirty() const {
    execute_for_each_sub_cell([](const sub_cell &cell) {
        cell.set_dirty();
    });
}

void
CellContainer::execute_for_each_sub_cell(const std::function<void(const CellContainer::sub_cell &)> &function) const {
    const auto grainSize = sub_cells().size() / config().nThreads();
    auto worker = [&](std::size_t, CellContainer::sub_cells_t::const_iterator begin,
                      CellContainer::sub_cells_t::const_iterator end) {
        for (auto it = begin; it != end; ++it) {
            function(*it);
        }
    };

    const auto& executor = *config().executor();
    std::vector<std::function<void(std::size_t)>> executables;
    executables.reserve(config().nThreads());

    auto it = sub_cells().begin();
    for (auto i = 0; i < config().nThreads() - 1; ++i) {
        executables.push_back(executor.pack(worker, it, it+grainSize));
        it += grainSize;
    }
    executables.push_back(executor.pack(worker, it, sub_cells().end()));
    executor.execute_and_wait(std::move(executables));

}

void CellContainer::update_root_size() {
    _root_size = vec3(_context.getBoxSize());
    _size = vec3(_context.getBoxSize());
}

CellContainer::~CellContainer() = default;

}
}
}
}
