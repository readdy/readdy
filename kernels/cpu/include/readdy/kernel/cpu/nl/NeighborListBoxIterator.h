/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * @file NeighborListBoxIterator.h
 * @brief << brief description >>
 * @author clonker
 * @date 12/11/17
 */

#pragma once

#include "CellLinkedList.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {



}
}
}
}

/*const auto &head = ccll.head();
const auto &list = ccll.list();

_elements.resize(1);
auto &neighbors = _elements.at(0);
neighbors.clear();

for (std::size_t cellIndex = 0; cellIndex < cix.size(); ++cellIndex) {

auto pptr = (*head.at(cellIndex)).load();
while (pptr != 0) {
auto pidx = pptr - 1;
const auto &entry = data.entry_at(pidx);
if(!entry.deactivated) {
neighbors.push_back(pidx);
auto n_neighbors_index = neighbors.size();
neighbors.push_back(0_z);
auto n_neighbors = 0_z;

{
auto ppptr = (*head.at(cellIndex)).load();
while (ppptr != 0) {
auto ppidx = ppptr - 1;
if (ppidx != pidx) {
const auto &pp = data.entry_at(ppidx);
if (!pp.deactivated) {
const auto distSquared = d2(entry.pos, pp.pos);
if (distSquared < cutoffSquared) {
neighbors.push_back(ppidx);
++n_neighbors;
}
}
}
ppptr = list.at(ppptr);
}
}

for (auto itNeighborCell = ccll.neighborsBegin(cellIndex);
itNeighborCell != ccll.neighborsEnd(cellIndex); ++itNeighborCell) {

auto nptr = (*head.at(*itNeighborCell)).load();
while (nptr != 0) {
auto nidx = nptr - 1;

const auto &neighbor = data.entry_at(nidx);
if (!neighbor.deactivated) {
const auto distSquared = d2(entry.pos, neighbor.pos);
if (distSquared < cutoffSquared) {
neighbors.push_back(nidx);
++n_neighbors;
}
}

nptr = list.at(nptr);
}
}
if (n_neighbors > 0) neighbors.at(n_neighbors_index) = n_neighbors;
}
pptr = list.at(pptr);
}
}*/