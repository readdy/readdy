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
 * @file CellLinkedList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu/nl/CellLinkedList.h>
#include <readdy/common/numeric.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

CellLinkedList::CellLinkedList(data_type &data, const readdy::model::Context &context,
                               const readdy::util::thread::Config &config)
        : _data(data), _context(context), _config(config) {}


void CellLinkedList::setUp(scalar skin, std::uint8_t radius, const util::PerformanceNode &node) {
    if(!_is_set_up || _skin != skin || _radius != radius) {
        auto t = node.timeit();

        _skin = skin;
        _radius = radius;
        _max_cutoff = _context.get().calculateMaxCutoff();
        if(_max_cutoff > 0) {
            auto size = _context.get().boxSize();
            auto desiredWidth = static_cast<scalar>((_max_cutoff + _skin)/static_cast<scalar>(radius));
            std::array<std::size_t, 3> dims{};
            for (int i = 0; i < 3; ++i) {
                dims[i] = static_cast<unsigned int>(std::max(1., std::floor(size[i] / desiredWidth)));
                _cellSize[i] = size[i] / static_cast<scalar>(dims[i]);
            }

            _cellIndex = util::Index3D(dims[0], dims[1], dims[2]);

            {
                // set up cell adjacency list
                auto t2 = node.subnode("setUpCellNeighbors").timeit();
                std::array<std::size_t, 3> nNeighbors{{_cellIndex[0], _cellIndex[1], _cellIndex[2]}};
                for (int i = 0; i < 3; ++i) {
                    nNeighbors[i] = std::min(nNeighbors[i], static_cast<std::size_t>(2 * radius + 1));
                }
                auto nAdjacentCells = nNeighbors[0] * nNeighbors[1] * nNeighbors[2];
                _cellNeighbors = util::Index2D(_cellIndex.size(), 1 + nAdjacentCells);
                _cellNeighborsContent.resize(_cellNeighbors.size());
                {
                    auto pbc = _context.get().periodicBoundaryConditions();
                    int r = _radius;
                    // local adjacency
                    std::vector<std::size_t> adj;
                    adj.reserve(1 + nAdjacentCells);
                    for (int i = 0; i < _cellIndex[0]; ++i) {
                        for (int j = 0; j < _cellIndex[1]; ++j) {
                            for (int k = 0; k < _cellIndex[2]; ++k) {
                                auto cellIdx = _cellIndex(static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                                                          static_cast<std::size_t>(k));

                                adj.resize(0);

                                for (int ii = i - r; ii <= i + r; ++ii) {
                                    for (int jj = j - r; jj <= j + r; ++jj) {
                                        for (int kk = k - r; kk <= k + r; ++kk) {
                                            if (ii == i && jj == j && kk == k) continue;
                                            auto adj_x = ii;
                                            if(pbc[0]) {
                                                auto cix = static_cast<int>(_cellIndex[0]);
                                                adj_x = (adj_x % cix + cix) % cix;
                                            }
                                            auto adj_y = jj;
                                            if(pbc[1]) {
                                                auto cix = static_cast<int>(_cellIndex[1]);
                                                adj_y = (adj_y % cix + cix) % cix;
                                            }
                                            auto adj_z = kk;
                                            if(pbc[2]) {
                                                auto cix = static_cast<int>(_cellIndex[2]);
                                                adj_z = (adj_z % cix + cix) % cix;
                                            }

                                            if (adj_x >= 0 && adj_y >= 0 && adj_z >= 0
                                                && adj_x < _cellIndex[0] && adj_y < _cellIndex[1] &&
                                                adj_z < _cellIndex[2]) {
                                                adj.push_back(_cellIndex(adj_x, adj_y, adj_z));
                                            }
                                        }
                                    }
                                }

                                std::sort(adj.begin(), adj.end());
                                adj.erase(std::unique(std::begin(adj), std::end(adj)), std::end(adj));

                                auto begin = _cellNeighbors(cellIdx, 0_z);
                                _cellNeighborsContent[begin] = adj.size();
                                std::copy(adj.begin(), adj.end(), &_cellNeighborsContent.at(begin + 1));
                            }
                        }
                    }
                }
            }

            if (_max_cutoff > 0) {
                setUpBins(node.subnode("setUpBins"));
            }

            _is_set_up = true;
        }
    }
}

const util::Index3D &CellLinkedList::cellIndex() const {
    return _cellIndex;
}

const util::Index2D &CellLinkedList::neighborIndex() const {
    return _cellNeighbors;
}

const std::vector<std::size_t> &CellLinkedList::neighbors() const {
    return _cellNeighborsContent;
}

std::size_t *CellLinkedList::neighborsBegin(std::size_t cellIndex) {
    auto beginIdx = _cellNeighbors(cellIndex, 1_z);
    return &_cellNeighborsContent.at(beginIdx);
}

const std::size_t *CellLinkedList::neighborsBegin(std::size_t cellIndex) const {
    auto beginIdx = _cellNeighbors(cellIndex, 1_z);
    return &_cellNeighborsContent.at(beginIdx);
}

std::size_t *CellLinkedList::neighborsEnd(std::size_t cellIndex) {
    return neighborsBegin(cellIndex) + nNeighbors(cellIndex);
}

const std::size_t *CellLinkedList::neighborsEnd(std::size_t cellIndex) const {
    return neighborsBegin(cellIndex) + nNeighbors(cellIndex);
}

std::size_t CellLinkedList::nNeighbors(std::size_t cellIndex) const {
    return _cellNeighborsContent.at(_cellNeighbors(cellIndex, 0_z));
}

CellLinkedList::data_type &CellLinkedList::data() {
    return _data.get();
}

const CellLinkedList::data_type &CellLinkedList::data() const {
    return _data.get();
}

scalar CellLinkedList::maxCutoff() const {
    return _max_cutoff;
}

CompactCellLinkedList::CompactCellLinkedList(data_type &data, const readdy::model::Context &context,
                                             const util::thread::Config &config) : CellLinkedList(data, context,
                                                                                                  config) {}

void CompactCellLinkedList::update(const util::PerformanceNode &node) {
    auto t = node.timeit();
    setUpBins(node.subnode("setUpBins"));
}

void CompactCellLinkedList::clear() {
    _head.resize(0);
    _list.resize(0);
}

std::size_t *CompactCellLinkedList::particlesBegin(std::size_t /*cellIndex*/) {
    throw std::logic_error("no particlesBegin for CompactCLL");
}

const std::size_t *CompactCellLinkedList::particlesBegin(std::size_t /*cellIndex*/) const {
    throw std::logic_error("no particlesBegin for CompactCLL");
}

std::size_t *CompactCellLinkedList::particlesEnd(std::size_t /*cellIndex*/) {
    throw std::logic_error("no particlesEnd for CompactCLL");
}

const std::size_t *CompactCellLinkedList::particlesEnd(std::size_t /*cellIndex*/) const {
    throw std::logic_error("no particlesEnd for CompactCLL");
}

std::size_t CompactCellLinkedList::nParticles(std::size_t cellIndex) const {
    std::size_t size = 0;
    auto pidx = (*_head.at(cellIndex)).load();
    while(pidx != 0) {
        pidx = _list.at(pidx);
        ++size;
    }
    return size;
}

template<>
void CompactCellLinkedList::fillBins<true>(const util::PerformanceNode &node) {
    auto t = node.timeit();
    auto boxSize = _context.get().boxSize();
    std::size_t pidx = 1;
    for (const auto &entry : _data.get()) {
        if (!entry.deactivated) {
            const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
            const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
            const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
            const auto cellIndex = _cellIndex(i, j, k);
            _list[pidx] = *_head.at(cellIndex);
            *_head[cellIndex] = pidx;
        }
        ++pidx;
    }
}

template<>
void CompactCellLinkedList::fillBins<false>(const util::PerformanceNode &node) {
    auto t = node.timeit();
    auto boxSize = _context.get().boxSize();
    const auto &data = _data.get();
    const auto grainSize = data.size() / _config.get().nThreads();
    const auto &executor = *_config.get().executor();
    const auto cellSize = _cellSize;
    const auto &cellIndex = _cellIndex;

    auto &list = _list;
    auto &head = _head;

    auto worker = [&data, &cellIndex, cellSize, &list, &head, boxSize]
            (std::size_t tid, std::size_t begin_pidx, std::size_t end_pidx) {
        auto it = data.begin() + begin_pidx - 1;
        auto pidx = begin_pidx;
        while(it != data.begin() + end_pidx - 1) {
            const auto &entry = *it;
            if (!entry.deactivated) {
                const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / cellSize.x));
                const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / cellSize.y));
                const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / cellSize.z));
                const auto cix = cellIndex(i, j, k);
                auto &atomic = *head.at(cix);
                // perform CAS
                auto currentHead = atomic.load();
                while(!atomic.compare_exchange_weak(currentHead, pidx)) {}
                list[pidx] = currentHead;
            }
            ++pidx;
            ++it;
        }
    };

    std::vector<std::function<void(std::size_t)>> executables;
    executables.reserve(_config.get().nThreads());
    auto it = 1_z;
    for (int i = 0; i < _config.get().nThreads() - 1; ++i) {
        executables.push_back(executor.pack(worker, it, it + grainSize));
        it += grainSize;
    }
    executables.push_back(executor.pack(worker, it, data.size()+1));
    executor.execute_and_wait(std::move(executables));
}

void CompactCellLinkedList::setUpBins(const util::PerformanceNode &node) {
    if(_max_cutoff > 0) {
        auto t = node.timeit();
        {
            auto tt = node.subnode("allocate").timeit();
            auto nParticles = _data.get().size();
            _head.clear();
            _head.resize(_cellIndex.size());
            _list.resize(0);
            _list.resize(nParticles + 1);
        }
        if(_serial) {
            fillBins<true>(node.subnode("fillBins serial"));
        } else {
            fillBins<false>(node.subnode("fillBins parallel"));
        }
    }
}

const CompactCellLinkedList::HEAD &CompactCellLinkedList::head() const {
    return _head;
}

const std::vector<std::size_t> &CompactCellLinkedList::list() const {
    return _list;
}

bool &CompactCellLinkedList::serial() {
    return _serial;
}

const bool &CompactCellLinkedList::serial() const {
    return _serial;
}

void CompactCellLinkedList::forEachParticlePair(const pair_callback &f) const {
    const auto &cix = cellIndex();
    const auto &data = _data.get();

    for (std::size_t cellIndex = 0; cellIndex < cix.size(); ++cellIndex) {

        auto pptr = (*_head.at(cellIndex)).load();
        while (pptr != 0) {
            auto pidx = pptr - 1;
            const auto &entry = data.entry_at(pidx);
            if(!entry.deactivated) {
                {
                    auto ppptr = (*_head.at(cellIndex)).load();
                    while (ppptr != 0) {
                        auto ppidx = ppptr - 1;
                        if (ppidx != pidx) {
                            const auto &pp = data.entry_at(ppidx);
                            if (!pp.deactivated) {
                                f(entry, pp);
                            }
                        }
                        ppptr = _list.at(ppptr);
                    }
                }

                for (auto itNeighborCell = neighborsBegin(cellIndex);
                     itNeighborCell != neighborsEnd(cellIndex); ++itNeighborCell) {

                    auto nptr = (*_head.at(*itNeighborCell)).load();
                    while (nptr != 0) {
                        auto nidx = nptr - 1;

                        const auto &neighbor = data.entry_at(nidx);
                        if (!neighbor.deactivated) {
                            f(entry, neighbor);
                        }

                        nptr = _list.at(nptr);
                    }
                }
            }
            pptr = _list.at(pptr);
        }
    }
}

void CompactCellLinkedList::forEachParticlePairParallel(const pair_callback &f) const {
    const auto &cix = cellIndex();
    const auto &data = _data.get();


    const auto grainSize = cix.size() / _config.get().nThreads();
    auto worker = [this, cix, &data, &f](std::size_t tid, std::size_t begin, std::size_t end) {
        const auto &head = this->head();
        const auto &list = this->list();
        for (std::size_t cellIndex = begin; cellIndex < end; ++cellIndex) {

            auto pptr = (*head.at(cellIndex)).load();
            while (pptr != 0) {
                auto pidx = pptr - 1;
                const auto &entry = data.entry_at(pidx);
                if(!entry.deactivated) {
                    {
                        auto ppptr = (*head.at(cellIndex)).load();
                        while (ppptr != 0) {
                            auto ppidx = ppptr - 1;
                            if (ppidx != pidx) {
                                const auto &pp = data.entry_at(ppidx);
                                if (!pp.deactivated) {
                                    f(entry, pp);
                                }
                            }
                            ppptr = list.at(ppptr);
                        }
                    }

                    for (auto itNeighborCell = neighborsBegin(cellIndex);
                         itNeighborCell != neighborsEnd(cellIndex); ++itNeighborCell) {

                        auto nptr = (*head.at(*itNeighborCell)).load();
                        while (nptr != 0) {
                            auto nidx = nptr - 1;

                            const auto &neighbor = data.entry_at(nidx);
                            if (!neighbor.deactivated) {
                                f(entry, neighbor);
                            }

                            nptr = list.at(nptr);
                        }
                    }
                }
                pptr = list.at(pptr);
            }
        }
    };
    {
        const auto &executor = *_config.get().executor();
        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(_config.get().nThreads());
        auto it = 0_z;
        for (int i = 0; i < _config.get().nThreads() - 1; ++i) {
            executables.push_back(executor.pack(worker, it, it + grainSize));
            it += grainSize;
        }
        executables.push_back(executor.pack(worker, it, cix.size()));
        executor.execute_and_wait(std::move(executables));
    }
}

BoxIterator CompactCellLinkedList::cellParticlesBegin(std::size_t cellIndex) const {
    auto head = (*_head.at(cellIndex)).load();
    return {*this, head};
}

BoxIterator CompactCellLinkedList::cellParticlesEnd(std::size_t) const {
    return {*this, 0};
}

NeighborsIterator CompactCellLinkedList::cellNeighborsBegin(std::size_t cellIndex) const {
    auto ix = (*_head.at(cellIndex)).load();
    if(ix > 0) {
        std::size_t nix = _list.at(ix);
        return {*this, cellIndex, nullptr, ix, nix};
    } else {
        return {*this, cellIndex, nullptr, 0, 0};
    }
}

NeighborsIterator CompactCellLinkedList::cellNeighborsEnd(std::size_t cellIndex) const {
    return {*this, cellIndex, neighborsEnd(cellIndex), 0, 0};
}

std::size_t CompactCellLinkedList::nCells() const {
    return _head.size();
}

bool CompactCellLinkedList::cellEmpty(std::size_t index) const {
    return (*_head.at(index)).load() == 0;
}

MacroBoxIterator CompactCellLinkedList::macroCellParticlesBegin(std::size_t cellIndex) const {
    return {*this, cellIndex, nullptr};
}

MacroBoxIterator CompactCellLinkedList::macroCellParticlesEnd(std::size_t cellIndex) const {
    return {*this, cellIndex, nullptr, true};
}

BoxIterator::BoxIterator(const CompactCellLinkedList &ccll, std::size_t state)
        : _ccll(ccll), _state(state) { }

BoxIterator &BoxIterator::operator++() {
    if (_state != 0) {
        _state = _ccll.get().list().at(_state);
    }
    return *this;
}

BoxIterator::value_type BoxIterator::operator*() const  {
    return _state - 1;
}

bool BoxIterator::operator==(const BoxIterator &rhs) const {
    return _state == rhs._state;
}

bool BoxIterator::operator!=(const BoxIterator &rhs) const {
    return !(*this == rhs);
}

MacroBoxIterator::MacroBoxIterator(const CompactCellLinkedList &ccll, std::size_t centerCell, 
                                   const std::size_t *currentCell, bool end)
        : _ccll(ccll), _centerCell(centerCell),
          _neighborCellsBegin(ccll.neighborsBegin(centerCell)), _neighborCellsEnd(ccll.neighborsEnd(centerCell)),
          _currentCell(getCurrentCell(currentCell)),
          _currentBoxIt(getCurrentBoxIterator()), _currentBoxEnd(ccll, 0) {
    
    if(_currentCell != _neighborCellsEnd && !end) {
        _state = *_currentBoxIt;
    } else {
        _state = 0;
        _currentBoxIt = _currentBoxEnd;
    }
}

MacroBoxIterator &MacroBoxIterator::operator++() {
    ++_currentBoxIt;
    if(_currentBoxIt == _currentBoxEnd) {
        if(_currentCell != _neighborCellsEnd) {
            if(*_currentCell == _centerCell) {
                _currentCell = getCurrentCell(_neighborCellsBegin);
            } else {
                _currentCell = getCurrentCell(_currentCell+1);
            }
        }
        _currentBoxIt = getCurrentBoxIterator();
    }
    if(_currentBoxIt != _currentBoxEnd) {
        _state = *_currentBoxIt;
    }
    return *this;
}

const std::size_t *MacroBoxIterator::getCurrentCell(const std::size_t *initial) const {
    if(initial == nullptr) initial = &_centerCell;
    const auto &ccll = _ccll.get();
    auto result = initial;
    if(ccll.cellEmpty(*result)) {
        if(*result == _centerCell) {
            result = _neighborCellsBegin;
        } else if (result != _neighborCellsEnd) {
            result = initial+1;
        }
        while(result != _neighborCellsEnd && ccll.cellEmpty(*result)) {
            ++result;
        }
    }
    return result;
}

BoxIterator MacroBoxIterator::getCurrentBoxIterator() const {
    const auto &ccll = _ccll.get();
    if(_currentCell != _neighborCellsEnd) {
        const auto head = (*ccll.head().at(*_currentCell)).load();
        return {ccll, head};
    } else {
        return {ccll, 0};
    }
}

const size_t &MacroBoxIterator::operator*() const {
    return _state;
}

MacroBoxIterator::pointer MacroBoxIterator::operator->() const {
    return &_state;
}

bool MacroBoxIterator::operator==(const MacroBoxIterator &rhs) const {
    return _currentBoxIt == rhs._currentBoxIt;
}

bool MacroBoxIterator::operator!=(const MacroBoxIterator &rhs) const {
    return _currentBoxIt != rhs._currentBoxIt;
}


NeighborsIterator::NeighborsIterator(const CompactCellLinkedList &ccll, std::size_t cell, const std::size_t *cellAt,
                                     std::size_t state, std::size_t stateNeigh)
        : _ccll(ccll), _state(state), _stateNeigh(stateNeigh), _cellAt(cellAt), _cell(cell),
          _innerState(ccll, stateNeigh), _outerState(ccll, state), _innerStateEnd(ccll, 0),
          _outerStateEnd(ccll.cellParticlesEnd(cell)) {
    _neighborCellsBegin = _ccll.get().neighborsBegin(cell);
    _neighborCellsEnd = _ccll.get().neighborsEnd(cell);

    // this deals with the special case of no neighboring cells
    if(_neighborCellsBegin == _neighborCellsEnd) {
        _neighborCellsBegin = nullptr;
        _neighborCellsEnd = nullptr;
    }

    if(_innerState != _innerStateEnd && _outerState != _outerStateEnd) {
        _currentValue = std::make_tuple(*_outerState, *_innerState);
    } else {
        if(_outerState == _outerStateEnd) {
            _innerState = _innerStateEnd;
            _outerState = _outerStateEnd;
        } else {
            // inner state is at end
            if(_cellAt == nullptr) {
                _cellAt = _neighborCellsBegin;
            }
            if(_cellAt != nullptr) {
                ++_cellAt;
            }
            while(_cellAt != _neighborCellsEnd && _ccll.get().cellEmpty(*_cellAt)) {
                ++_cellAt;
            }
            if(_cellAt != _neighborCellsEnd) {
                _innerState = _ccll.get().cellParticlesBegin(*_cellAt);
            } else {
                _innerState = _innerStateEnd;
                _outerState = _outerStateEnd;
            }
        }
    }
}

NeighborsIterator &NeighborsIterator::operator++() {
    // the outer loop through all particles of the inner cell
    if(_cellAt == _neighborCellsEnd && _innerState == _innerStateEnd) {
        // the last cell has been reached, begin anew with the next particle in the own cell
        ++_outerState;
        _cellAt = nullptr;
        if(_outerState != _outerStateEnd) {
            _innerState = _ccll.get().cellParticlesBegin(_cell);
        } else {
            // once the end is reached, simply set inner state to end
            _innerState = _innerStateEnd;
        }
    } else if(_cellAt == nullptr) {
        // we are inside "own" cell, so check if we have ourself as particle in the inner state or need to advance
        // the neighboring cell
        ++_innerState;
        if(_innerState != _innerStateEnd && _innerState == _outerState) {
            // inner and outer point to the same particle, skip
            ++_innerState;
        }

        if(_innerState == _innerStateEnd) {
            // the end is reached, proceed (if available) to next cell
            if(_neighborCellsBegin != nullptr) {
                _cellAt = _neighborCellsBegin;
                while(_cellAt != _neighborCellsEnd && _ccll.get().cellEmpty(*_cellAt)) {
                    ++_cellAt;
                }
                if(_cellAt != _neighborCellsEnd) {
                    _innerState = _ccll.get().cellParticlesBegin(*_cellAt);
                } else {
                    ++_outerState;
                    if(_outerState != _outerStateEnd) {
                        _cellAt = nullptr;
                        _innerState = _ccll.get().cellParticlesBegin(_cell);
                    } else {
                        _innerState = _innerStateEnd;
                    }
                }
            }
        }
    } else {
        // we are in one of the surrounding cells so the particle pointed to by outerState cannot be in here
        ++_innerState;

        if(_innerState == _innerStateEnd) {
            // proceed to next cell
            ++_cellAt;
            while(_cellAt != _neighborCellsEnd && _ccll.get().cellEmpty(*_cellAt)) {
                ++_cellAt;
            }
            if(_cellAt != _neighborCellsEnd) {
                _innerState = _ccll.get().cellParticlesBegin(*_cellAt);
                _innerStateEnd = _ccll.get().cellParticlesEnd(*_cellAt);
            } else {
                _innerState = {_ccll.get(), 0};
            }
        }
    }
    if(_innerState != _innerStateEnd && _outerState != _outerStateEnd) {
        _currentValue = std::make_tuple(*_outerState, *_innerState);
    }
    return *this;
}

NeighborsIterator::reference NeighborsIterator::operator*() const {
    return _currentValue;
}

bool NeighborsIterator::operator==(const NeighborsIterator &rhs) const {
    return _outerState == rhs._outerState && _innerState == rhs._innerState;
}

bool NeighborsIterator::operator!=(const NeighborsIterator &rhs) const {
    return !(*this == rhs);
}


}
}
}
}
