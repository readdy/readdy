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
 * @file CLLNeighborList.cpp
 * @brief 
 * @author clonker
 * @date 9/13/17
 */

#include <readdy/kernel/cpu_legacy/nl/CLLNeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {


ContiguousCLLNeighborList::ContiguousCLLNeighborList(std::uint8_t cll_radius,
                                                     const readdy::model::Context &context,
                                                     const readdy::util::thread::Config &config)
        : NeighborList(context, config), _data(context, config), ccll(_data, context, config),
          ccllContainer({context, config, _data}, ccll), cll_radius(cll_radius) {}


ContiguousCLLNeighborList::ContiguousCLLNeighborList(data::EntryDataContainer *data, std::uint8_t cll_radius,
                                                     const readdy::model::Context &context,
                                                     const readdy::util::thread::Config &config)
        : NeighborList(context, config), _data(data), ccll(_data, context, config),
          ccllContainer({context, config, _data}, ccll), cll_radius(cll_radius) {}

void ContiguousCLLNeighborList::set_up(const util::PerformanceNode &node) {
    auto t = node.timeit();
    _max_cutoff = _context.get().calculateMaxCutoff();
    _max_cutoff_skin_squared = (_max_cutoff + _skin) * (_max_cutoff + _skin);

    ccll.setUp(_skin, cll_radius, node.subnode("setUp CLL"));
    fill_verlet_list(node.subnode("fill verlet list"));
    _is_set_up = true;
}

void ContiguousCLLNeighborList::fill_verlet_list(const util::PerformanceNode &node) {
    if (_max_cutoff > 0) {
        ccllContainer.update(_max_cutoff_skin_squared, node.subnode("update neighbor list container"));
    }
}

void ContiguousCLLNeighborList::update(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (!_is_set_up) {
        set_up(node.subnode("setUp"));
    } else {
        ccll.update(node.subnode("update CLL"));
        fill_verlet_list(node.subnode("fill verlet list"));
    }
}

void ContiguousCLLNeighborList::clear(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (_max_cutoff > 0) {
        ccll.clear();
        ccllContainer.clear();
    }
}

void ContiguousCLLNeighborList::updateData(DataUpdate &&update) {
    _data.update(std::forward<DataUpdate>(update));
}

bool ContiguousCLLNeighborList::is_adaptive() const {
    return false;
}

NeighborList::const_iterator ContiguousCLLNeighborList::cbegin() const {
    return NeighborListIterator{ccllContainer.begin(), ccllContainer.end(), false};
}

NeighborList::const_iterator ContiguousCLLNeighborList::cend() const {
    return NeighborListIterator{ccllContainer.end(), ccllContainer.end(), false};
}

const data::EntryDataContainer *ContiguousCLLNeighborList::data() const {
    return &_data;
}

data::EntryDataContainer *ContiguousCLLNeighborList::data() {
    return &_data;
}

std::size_t ContiguousCLLNeighborList::size() const {
    return ccllContainer.size();
}

DynamicCLLNeighborList::DynamicCLLNeighborList(data::EntryDataContainer *data,
                                               const readdy::model::Context &context,
                                               const readdy::util::thread::Config &config)
        : NeighborList(context, config), _data(data), dcll(_data, context, config),
          dcllContainer({context, config, _data}, dcll) {}


DynamicCLLNeighborList::DynamicCLLNeighborList(const readdy::model::Context &context,
                                               const readdy::util::thread::Config &config)
        : NeighborList(context, config), _data(context, config), dcll(_data, context, config),
          dcllContainer({context, config, _data}, dcll) {}

void DynamicCLLNeighborList::set_up(const util::PerformanceNode &node) {
    auto t = node.timeit();
    _max_cutoff = _context.get().calculateMaxCutoff();
    _max_cutoff_skin_squared = (_max_cutoff + _skin) * (_max_cutoff + _skin);

    dcll.setUp(_skin, 1, node.subnode("setUp CLL"));
    fill_verlet_list(node.subnode("fill verlet list"));
    _is_set_up = true;
}

void DynamicCLLNeighborList::update(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (!_is_set_up) {
        set_up(node.subnode("setUp"));
    } else {
        dcll.update(node.subnode("update CLL"));
        fill_verlet_list(node.subnode("fill verlet list"));
    }
}

void DynamicCLLNeighborList::clear(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (_max_cutoff > 0) {
        dcll.clear();
        dcllContainer.clear();
    }
}

void DynamicCLLNeighborList::updateData(DataUpdate &&update) {
    _data.update(std::forward<DataUpdate>(update));
}

void DynamicCLLNeighborList::fill_verlet_list(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (_max_cutoff > 0) {
        dcllContainer.update(_max_cutoff_skin_squared, node.subnode("fill verlet list"));
    }
}

bool DynamicCLLNeighborList::is_adaptive() const {
    return false;
}

NeighborList::const_iterator DynamicCLLNeighborList::cbegin() const {
    return NeighborListIterator{dcllContainer.begin(), dcllContainer.end(), false};
}

NeighborList::const_iterator DynamicCLLNeighborList::cend() const {
    return NeighborListIterator{dcllContainer.end(), dcllContainer.end(), false};
}

const data::EntryDataContainer *DynamicCLLNeighborList::data() const {
    return &_data;
}

data::EntryDataContainer *DynamicCLLNeighborList::data() {
    return &_data;
}

std::size_t DynamicCLLNeighborList::size() const {
    return dcllContainer.size();
}

CompactCLLNeighborList::CompactCLLNeighborList(std::uint8_t cll_radius, const readdy::model::Context &context,
                                               const util::thread::Config &config)
        : NeighborList(context, config), _data(context, config), ccll(_data, context, config), cll_radius(cll_radius),
          ccllContainer({context, config, _data}, ccll) {}

CompactCLLNeighborList::CompactCLLNeighborList(data::EntryDataContainer *data, std::uint8_t cll_radius,
                                               const readdy::model::Context &context,
                                               const util::thread::Config &config)
        : NeighborList(context, config), _data(data), ccll(_data, context, config), cll_radius(cll_radius),
          ccllContainer({context, config, _data}, ccll) {}

void CompactCLLNeighborList::set_up(const util::PerformanceNode &node) {
    auto t = node.timeit();
    _max_cutoff = _context.get().calculateMaxCutoff();
    _max_cutoff_skin_squared = (_max_cutoff + _skin) * (_max_cutoff + _skin);

    ccll.setUp(_skin, cll_radius, node.subnode("setUp CLL"));
    fill_verlet_list(node.subnode("fill verlet list"));
    _is_set_up = true;
}

void CompactCLLNeighborList::update(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (!_is_set_up) {
        set_up(node.subnode("setUp"));
    } else {
        ccll.update(node.subnode("update CLL"));
        fill_verlet_list(node.subnode("fill verlet list"));
    }
}

void CompactCLLNeighborList::clear(const util::PerformanceNode &node) {
    ccll.clear();
}

void CompactCLLNeighborList::updateData(DataUpdate &&update) {
    _data.update(std::forward<DataUpdate>(update));
}

void CompactCLLNeighborList::fill_verlet_list(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (_max_cutoff > 0) {
        ccllContainer.update(_max_cutoff_skin_squared, node.subnode("fill verlet list"));
    }
}

bool CompactCLLNeighborList::is_adaptive() const {
    return false;
}

NeighborList::const_iterator CompactCLLNeighborList::cbegin() const {
    return NeighborListIterator{ccllContainer.begin(), ccllContainer.end(), false};
}

NeighborList::const_iterator CompactCLLNeighborList::cend() const {
    return NeighborListIterator{ccllContainer.end(), ccllContainer.end(), false};
}

const data::EntryDataContainer *CompactCLLNeighborList::data() const {
    return &_data;
}

data::EntryDataContainer *CompactCLLNeighborList::data() {
    return &_data;
}

const CompactCLLNeighborListContainer &CompactCLLNeighborList::container() const {
    return ccllContainer;
}

const CompactCellLinkedList &CompactCLLNeighborList::cellLinkedList() const {
    return ccll;
}

CompactCellLinkedList &CompactCLLNeighborList::cellLinkedList() {
    return ccll;
}

std::size_t CompactCLLNeighborList::size() const {
    return ccllContainer.size();
}

}
}
}
}