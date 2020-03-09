/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/

/**
 * « detailed description »
 *
 * @file MPINeighborList.h
 * @brief « brief description »
 * @author chrisfroe
 * @date 14.06.19
 */

#pragma once

#include <readdy/kernel/mpi/model/MPIParticleData.h>
#include <readdy/kernel/mpi/model/MPIDomain.h>

namespace readdy::kernel::mpi::model {

class BoxIterator;

class CellLinkedList {
public:
    using Data = readdy::kernel::mpi::MPIDataContainer;
    using CellRadius = std::uint8_t;
    using HEAD = std::unordered_map<std::size_t,std::size_t>;
    using LIST = std::vector<std::size_t>;
    using CellNeighbors = std::unordered_map<std::size_t, std::vector<std::size_t>>;
    using EntryCref = const Data::EntryType &;
    using PairCallback = std::function<void(EntryCref, EntryCref)>;

    using IteratorBounds = std::tuple<std::size_t, std::size_t>;

    CellLinkedList(Data &data, const readdy::model::Context &context)
            : _data(data), _context(context), _head{}, _list{}, _radius{0} {};

    void setUp(const model::MPIDomain* domain) {
        if (!_isSetUp) {
            _domain = domain;

            // Each domain is subdivided into `cellsExtent[coord]` cells along axis coord
            // This guarantees that domain boundaries are also cell boundaries
            const auto boxSize = _context.get().boxSize();
            std::array<std::size_t, 3> nCellsPerAxis{};
            std::array<std::size_t, 3> cellsExtent{}; // number of cells per domain per axis
            std::array<std::size_t, 3> cellsOrigin{}; // ijk of this domain's origin cell, i.e. the lower left cell
            for (std::size_t coord = 0; coord < 3; ++coord) {
                cellsExtent[coord] = static_cast<std::size_t>(std::max(1., std::floor(domain->extent()[coord] / domain->haloThickness())));
                cellsOrigin[coord] = domain->myIdx()[coord] * cellsExtent[coord];
                nCellsPerAxis[coord] = cellsExtent[coord] * domain->nDomainsPerAxis()[coord];
                _cellSize[coord] = boxSize[coord] / static_cast<scalar>(nCellsPerAxis[coord]);
            }

            _cellIndex = readdy::util::Index3D(nCellsPerAxis[0], nCellsPerAxis[1], nCellsPerAxis[2]);

            {
                // todo
                // data structure that takes a map<uint,vector<uint>> as ctor arg, and provides a dense array access
                // via Index2D and a backing structure which is one vector (with padding)
                // difference to old impl here would be that the backing vector is sparse wrt to the index uint
                // i.e. internally there is a remapping umap<uint,uint> from sparse to dense indices
                // can this be efficient? due to umap memory is again fragmented
                // -> my hunch is that one would just need a lookup-performance-optimized unordered_map
                {
                    // local adjacency, iterate over cells in domain core
                    for (int i = cellsOrigin[0]; i < cellsOrigin[0] + cellsExtent[0]; ++i) {
                    for (int j = cellsOrigin[1]; j < cellsOrigin[1] + cellsExtent[1]; ++j) {
                    for (int k = cellsOrigin[2]; k < cellsOrigin[2] + cellsExtent[2]; ++k) {
                        const auto cellIdx = _cellIndex(i,j,k);
                        _cellsInCore.push_back(cellIdx);

                        // for di,dj,dk, this also reaches neighbor cells that overlap with halo
                        for (int di=-1; di<2; ++di) {
                        for (int dj=-1; dj<2; ++dj) {
                        for (int dk=-1; dk<2; ++dk) {
                            if (di==0 and dj==0 and dk==0) {
                                continue; // skip self
                            }
                            std::array<int, 3> otherCell = {i+di, j+dj, k+dk};

                            // fix boundaries
                            bool isValidCell = true;
                            for (std::uint8_t axis = 0; axis < 3; ++axis) {
                                auto nCells = static_cast<int>(_cellIndex[axis]);
                                const auto pbc = _context.get().periodicBoundaryConditions();
                                if (pbc[axis] && nCells > 2) {
                                    if (-1 <= otherCell[axis] and otherCell[axis] <= nCells) {
                                        otherCell.at(axis) = (otherCell.at(axis) % nCells + nCells) % nCells;
                                    } else {
                                        isValidCell = false;
                                    }
                                } else if (0 <= otherCell[axis] and otherCell[axis] < nCells) {
                                    // all good, cell is within boxSize
                                } else {
                                    isValidCell = false;
                                }
                            }

                            // add otherCell to neighborhood of cell, avoid double neighborliness
                            if (isValidCell) {
                                assert(std::all_of(otherCell.begin(), otherCell.end(), [](const auto& x){return x>=0;}));
                                const auto otherIdx = _cellIndex(otherCell[0],otherCell[1],otherCell[2]);
                                const auto cellCenter = Vec3(
                                        otherCell[0] * _cellSize[0] + 0.5 * _cellSize[0],
                                        otherCell[1] * _cellSize[1] + 0.5 * _cellSize[1],
                                        otherCell[2] * _cellSize[2] + 0.5 * _cellSize[2]);
                                if (domain->isInDomainCore(cellCenter)) {
                                    if (cellIdx < otherIdx) { // avoid double neighborliness for core cells
                                        _cellNeighbors.at(cellIdx).push_back(otherIdx);
                                    }
                                } else {
                                    // cells that are outside the core are only seen "from one side",
                                    // because we iterate over the core cells only.
                                    // thus always add them to the neighborhood
                                    _cellNeighbors.at(cellIdx).push_back(otherIdx);
                                    // additionally keep track of all cells that overlap with halo
                                    _cellsInHalo.push_back(otherIdx);
                                }
                            }

                        }
                        }
                        }
                    }
                    }
                    }
                }

                // clean up, why not
                std::sort(_cellsInCore.begin(), _cellsInCore.end());

                // remove duplicates out of _cellsInHalo
                std::sort(_cellsInHalo.begin(), _cellsInHalo.end());
                auto last = std::unique(std::begin(_cellsInHalo), std::end(_cellsInHalo));
                _cellsInHalo.erase(last, std::end(_cellsInHalo));
            }
            _isSetUp = true;
            setUpBins();
        }
    };

    void update() {
        setUpBins();
    }

    virtual void clear() {
        _head.clear();
        _list.resize(0);
        _isSetUp = false;
    }

    const readdy::util::Index3D &cellIndex() const {
        return _cellIndex;
    }

    const CellNeighbors &neighborIndex() const {
        return _cellNeighbors;
    }

    CellNeighbors::mapped_type::iterator neighborsBegin(std::size_t cellIndex) {
        return _cellNeighbors.at(cellIndex).begin();
    }

    CellNeighbors::mapped_type::const_iterator neighborsBegin(std::size_t cellIndex) const {
        return _cellNeighbors.at(cellIndex).cbegin();
    }

    CellNeighbors::mapped_type::iterator neighborsEnd(std::size_t cellIndex) {
        return _cellNeighbors.at(cellIndex).end();
    }

    CellNeighbors::mapped_type::const_iterator neighborsEnd(std::size_t cellIndex) const {
        return _cellNeighbors.at(cellIndex).cend();
    }

    [[nodiscard]] std::size_t nNeighbors(std::size_t cellIndex) const {
        return _cellNeighbors.at(cellIndex).size();
    }

    Data &data() {
        return _data.get();
    }

    const Data &data() const {
        return _data.get();
    }

    scalar cutoff() const {
        return _cutoff;
    }

    const HEAD &head() const {
        return _head;
    }

    const LIST &list() const {
        return _list;
    }

    const std::vector<std::size_t> &cellsInCore() const {
        return _cellsInCore;
    }

    const std::vector<std::size_t> &cellsInHalo() const {
        return _cellsInHalo;
    }

    template<typename Function>
    void forEachNeighbor(BoxIterator particle, const Function &function) const;

    template<typename Function>
    void forEachNeighbor(BoxIterator particle, std::size_t cell, const Function &function) const;

    bool cellEmpty(std::size_t index) const {
        return _head.at(index) == 0;
    }

    std::size_t cellOfParticle(std::size_t index) const {
        const auto &entry = data().entry_at(index);
        if (entry.deactivated) {
            throw std::invalid_argument("requested deactivated entry");
        }
        const auto &boxSize = _context.get().boxSize();
        if(!(-.5*boxSize[0] <= entry.pos.x && .5*boxSize[0] > entry.pos.x
             && -.5*boxSize[1] <= entry.pos.y && .5*boxSize[1] > entry.pos.y
             && -.5*boxSize[2] <= entry.pos.z && .5*boxSize[2] > entry.pos.z)) {
            throw std::logic_error("CellLinkedList: requested neighbors of particle that was out of bounds.");
        }
        const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
        const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
        const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
        return _cellIndex(i, j, k);
    }

    std::size_t nCells() const {
        return _cellIndex.size();
    }

    BoxIterator particlesBegin(std::size_t cellIndex);

    BoxIterator particlesBegin(std::size_t cellIndex) const;

    BoxIterator particlesEnd(std::size_t cellIndex);

    BoxIterator particlesEnd(std::size_t cellIndex) const;

    const Vec3 &cellSize() const {
        return _cellSize;
    }

protected:
    virtual void setUpBins() {
        if (_isSetUp) {
            auto nParticles = _data.get().size();
            _head.clear(); // head structure will be built lazily upon filling bins
            _list.resize(0);
            _list.resize(nParticles + 1); // _list[0] is terminator for a sequence of particles
            fillBins();
        } else {
            throw std::logic_error("Attempting to fill neighborlist bins, but cell structure is not set up yet");
        }
    }

    void fillBins() {
        const auto &boxSize = _context.get().boxSize();

        const auto particleInBox = [&boxSize](const Vec3 &pos) {
            return -.5*boxSize[0] <= pos.x && .5*boxSize[0] > pos.x
                   && -.5*boxSize[1] <= pos.y && .5*boxSize[1] > pos.y
                   && -.5*boxSize[2] <= pos.z && .5*boxSize[2] > pos.z;
        };

        std::size_t pidx = 1; // the list structure is 1-indexed, because 0 terminates the particle group
        for (const auto &entry : _data.get()) {
            if (not entry.deactivated) {
                if (particleInBox(entry.pos)) {
                    const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
                    const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
                    const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
                    const auto cellIndex = _cellIndex(i, j, k);
                    // if head[cellIndex] does not exist here it is default (0) which is the terminator
                    _list[pidx] = _head[cellIndex];
                    _head[cellIndex] = pidx;
                } else {
                    readdy::log::warn("Particle not in box, will not be contained in the neighbor-list");
                }
            }
            ++pidx;
        }
    }

    // head maps from cell indices to the first particle of a group in the list structure
    HEAD _head;
    // Linear list of particles, 1-indexed
    // Used to build a string of particles (that are contained in a cell),
    // index i refers to a particle index (i-1) and j=list[i] contains the index of the
    // next particle, this string is terminated if j=0
    LIST _list;


    bool _isSetUp{false};

    scalar _cutoff{0};
    std::uint8_t _radius;

    Vec3 _cellSize{0, 0, 0};

    readdy::util::Index3D _cellIndex;

    // maps from cell index to neighbor cell indices, consider a dense structure again
    CellNeighbors _cellNeighbors;

    // keep track which cells are in the core of the domain and which cells overlap with the halo region
    std::vector<std::size_t> _cellsInCore;
    std::vector<std::size_t> _cellsInHalo;

    std::reference_wrapper<Data> _data;
    std::reference_wrapper<const readdy::model::Context> _context;
    const model::MPIDomain * _domain{nullptr};
};

class BoxIterator {

    using alloc = std::allocator<std::size_t>;

public:

    using difference_type = typename alloc::difference_type;
    using value_type = typename alloc::value_type;
    using reference = typename alloc::const_reference;
    using pointer = typename alloc::const_pointer;
    using iterator_category = std::forward_iterator_tag;
    using size_type = CellLinkedList::LIST::size_type;

    BoxIterator(const CellLinkedList &ccll, std::size_t state) : _ccll(ccll), _state(state), _val(state - 1) {};

    BoxIterator(const BoxIterator &) = default;

    BoxIterator &operator=(const BoxIterator &) = delete;

    BoxIterator(BoxIterator &&) = default;

    BoxIterator &operator=(BoxIterator &&) = delete;

    ~BoxIterator() = default;

    const BoxIterator operator++(int) {
        BoxIterator tmp(*this);
        operator++();
        return tmp;
    }

    pointer operator->() const {
        return &_val;
    }

    BoxIterator &operator++() {
        _state = _ccll.list().at(_state);
        _val = _state - 1;
        return *this;
    }

    value_type operator*() const {
        return _val;
    }

    bool operator==(const BoxIterator &rhs) const {
        return _state == rhs._state;
    }

    bool operator!=(const BoxIterator &rhs) const {
        return _state != rhs._state;
    }

private:
    const CellLinkedList &_ccll;
    std::size_t _state, _val;
};

inline BoxIterator CellLinkedList::particlesBegin(std::size_t cellIndex) {
    return {*this, _head.at(cellIndex)};
}

inline BoxIterator CellLinkedList::particlesBegin(std::size_t cellIndex) const {
    return {*this, _head.at(cellIndex)};
}

inline BoxIterator CellLinkedList::particlesEnd(std::size_t /*cellIndex*/) const {
    return {*this, 0};
}

inline BoxIterator CellLinkedList::particlesEnd(std::size_t /*cellIndex*/) {
    return {*this, 0};
}

template<typename Function>
inline void CellLinkedList::forEachNeighbor(BoxIterator particle, const Function &function) const {
    forEachNeighbor(particle, cellOfParticle(*particle), function);
}

template<typename Function>
inline void CellLinkedList::forEachNeighbor(BoxIterator particle, std::size_t cell,
                                            const Function &function) const {
    std::for_each(std::next(particle, 1), particlesEnd(cell), [&function, particle](auto x) {
        function(x);
    });
    for (auto itNeighCell = neighborsBegin(cell); itNeighCell != neighborsEnd(cell); ++itNeighCell) {
        std::for_each(particlesBegin(*itNeighCell), particlesEnd(*itNeighCell), function);
    }
}

}
