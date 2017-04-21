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

#include <readdy/common/numeric.h>
#include <readdy/kernel/cpu/nl/CellContainer.h>
#include <readdy/kernel/cpu/nl/SubCell.h>
#include <readdy/kernel/cpu/model/CPUParticleData.h>
#include <readdy/kernel/cpu/util/config.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

void execute_for_each_sub_cell(const std::function<void(CellContainer::sub_cell&)>& fun,
                               CellContainer& container) {
    const auto grainSize = container.cells().size() / container.config().nThreads();
    auto worker = [&](CellContainer::sub_cells::iterator begin, CellContainer::sub_cells::iterator end) {
        for(auto it = begin; it != end; ++it) {
            fun(*it);
        }
    };
    std::vector<threading_model> threads;
    threads.reserve(container.config().nThreads());
    auto it = container.cells().begin();
    for(auto i = 0; i < container.config().nThreads() - 1; ++i) {
        threads.emplace_back(worker, it, it + grainSize);
        it += grainSize;
    }
    threads.emplace_back(worker, it, container.cells().end());
}

CellContainer::sub_cells &CellContainer::cells() {
    return _cells;
}

const CellContainer::sub_cells &CellContainer::cells() const {
    return _cells;
}

CellContainer::dimension &CellContainer::n_cells() {
    return _n_cells;
}

const CellContainer::dimension &CellContainer::n_cells() const {
    return _n_cells;
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

CellContainer::CellContainer(const model::CPUParticleData &data, const readdy::model::KernelContext& context,
                             const readdy::util::thread::Config& config)
        : _data(data), _context(context), _config(config) {}

CellContainer::sub_cell *const CellContainer::cell_for_position(const CellContainer::vec3 &pos) {
    const cell_index i = static_cast<const cell_index>(floor((pos.x - _offset.x + .5 * _size.x) / _cell_size.x));
    const cell_index j = static_cast<const cell_index>(floor((pos.y - _offset.y + .5 * _size.y) / _cell_size.y));
    const cell_index k = static_cast<const cell_index>(floor((pos.z - _offset.z + .5 * _size.z) / _cell_size.z));
    return cell_for_grid_index(std::tie(i, j, k));
}

const CellContainer::sub_cell *const CellContainer::cell_for_position(const CellContainer::vec3 &pos) const {
    const cell_index i = static_cast<const cell_index>(floor((pos.x - _offset.x + .5 * _size.x) / _cell_size.x));
    const cell_index j = static_cast<const cell_index>(floor((pos.y - _offset.y + .5 * _size.y) / _cell_size.y));
    const cell_index k = static_cast<const cell_index>(floor((pos.z - _offset.z + .5 * _size.z) / _cell_size.z));
    return cell_for_grid_index(std::tie(i, j, k));
}

CellContainer::sub_cell *const CellContainer::cell_for_grid_index(const CellContainer::grid_index &index) {
    auto i = std::get<0>(index);
    auto j = std::get<1>(index);
    auto k = std::get<2>(index);
    const auto &periodic = _context.getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, _n_cells[0]);
    else if (i < 0 || i >= _n_cells[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, _n_cells[1]);
    else if (j < 0 || j >= _n_cells[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, _n_cells[2]);
    else if (k < 0 || k >= _n_cells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, _n_cells[1], _n_cells[2]);
    if (cix < _cells.size()) {
        return &_cells.at(static_cast<cell_index>(cix));
    } else {
        log::critical("CPUNeighborList::getCell(nonconst): Requested cell ({},{},{})={}, but there are "
                              "only {} cells.", i, j, k, cix, _cells.size());
        throw std::runtime_error("tried to get cell index that was too large");
    }
}

const CellContainer::sub_cell *const CellContainer::cell_for_grid_index(const CellContainer::grid_index &index) const {
    auto i = std::get<0>(index);
    auto j = std::get<1>(index);
    auto k = std::get<2>(index);
    const auto &periodic = _context.getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, _n_cells[0]);
    else if (i < 0 || i >= _n_cells[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, _n_cells[1]);
    else if (j < 0 || j >= _n_cells[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, _n_cells[2]);
    else if (k < 0 || k >= _n_cells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, _n_cells[1], _n_cells[2]);
    if (cix < _cells.size()) {
        return &_cells.at(static_cast<cell_index>(cix));
    } else {
        log::critical("CPUNeighborList::getCell(const): Requested cell ({},{},{})={}, but there are "
                              "only {} cells.", i, j, k, cix, _cells.size());
        throw std::runtime_error("tried to get cell index that was too large");
    }
}

void CellContainer::update_displacements() {
    update_sub_cell_displacements();
    _maximal_displacements[0] = 0;
    _maximal_displacements[1] = 0;
    for(auto& cell : _cells) {
        const auto& updated_displacements = cell.maximal_displacements();
        if(updated_displacements[0] > _maximal_displacements[0]) {
            _maximal_displacements[0] = updated_displacements[0];
            _maximal_displacements[1] = std::max(_maximal_displacements[1], updated_displacements[1]);
        } else if(updated_displacements[0] > _maximal_displacements[1]) {
            _maximal_displacements[1] = updated_displacements[0];
        }
    }
}

void CellContainer::subdivide(const scalar desired_cell_width) {
    for (std::uint8_t i = 0; i < 3; ++i) {
        _n_cells[i] = std::max(static_cast<cell_index>(1), static_cast<cell_index>(floor(_size[i] / desired_cell_width)));
        _cell_size[i] = _size[i] / static_cast<double>(_n_cells[i]);
    }
    log::debug("resulting cell size = {}", _cell_size);
    _cells.reserve(_n_cells[0] * _n_cells[1] * _n_cells[2]);
    for (cell_index i = 0; i < _n_cells[0]; ++i) {
        for (cell_index j = 0; j < _n_cells[1]; ++j) {
            for (cell_index k = 0; k < _n_cells[2]; ++k) {
                _cells.emplace_back(*this, vec3{i * _cell_size.x, j * _cell_size.y, k * _cell_size.z});
                _cells.back()._size = _cell_size;
                _cells.back()._contiguous_index = get_contiguous_index(i, j, k, _n_cells[1], _n_cells[2]);
            }
        }
    }
    // todo set up neighbors and this method should be hidden in subcell (too much flexibility that we dont need right now)
}

const model::CPUParticleData &CellContainer::data() const {
    return _data;
}

const readdy::model::KernelContext &CellContainer::context() const {
    return _context;
}

void CellContainer::update_sub_cell_displacements() {
    execute_for_each_sub_cell([](sub_cell& cell) {
        cell.update_displacements();
    }, *this);
}

const readdy::util::thread::Config &CellContainer::config() const {
    return _config;
}


}
}
}
}