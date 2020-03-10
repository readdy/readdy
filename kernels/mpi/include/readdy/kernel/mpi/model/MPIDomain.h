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
 * Each worker/rank is responsible for one region in space and the particles that live in it. The region of space
 * is defined by the ranks' MPIDomain. It also provides neighborship information like: ranks of the adjacent domains.
 * Although this holds much valuable information for communicating amongst processes,
 * MPIDomain makes no calls to the MPI library, which allows for intensive testing/debugging without the MPI context.
 *
 * @file MPIDomain.h
 * @brief Spatial setup of MPI domains
 * @author chrisfroe
 * @date 17.06.19
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/model/Context.h>

namespace readdy::kernel::mpi::model {

class MPIDomain {
    /** Members which the master rank also knows */
public:

    [[nodiscard]] int rank() const { return _rank; }
    [[nodiscard]] int worldSize() const { return _worldSize; }
    [[nodiscard]] scalar haloThickness() const { return _haloThickness; }

private:

    int _rank;
    int _worldSize;
    scalar _haloThickness;

    int _nUsedRanks; // counts master rank and all workers
    int _nWorkerRanks; // counts all workers, which is a subset of used ranks
    int _nIdleRanks; // counts all idle
    std::vector<int> _workerRanks; // can be returned as const ref to conveniently iterate over ranks
    bool _isNotIdle{true};
    std::array<std::size_t, 3> _nDomainsPerAxis{};
    readdy::util::Index3D _domainIndex; // rank of (ijk) is domainIndex(i,j,k)+1

    /** The following members will only be defined for rank != 0 */

    // origin and extent define the core region of the domain
    Vec3 _origin; // lower-left corner
    Vec3 _extent;
    // these define the extended region core+halo
    Vec3 _originWithHalo;
    Vec3 _extentWithHalo;

    std::array<std::size_t, 3> _myIdx{}; // (ijk) indices of this domain
public:
    // Neighborhood stuff
    // map from ([0-2], [0-2], [0-2]) to the index of the 27 neighbors (including self)
    const readdy::util::Index3D neighborIndex = readdy::util::Index3D(static_cast<std::size_t>(3), static_cast<std::size_t>(3),
                                                      static_cast<std::size_t>(3));

    enum NeighborType {
        self, // is a valid neighbor, but due to periodicity it is this domain
        nan, // not a neighbor (i.e. end of simulation box and no periodic bound.)
        regular // neighbor is another domain ()
    };

private:
    // the backing structures that the _neighborIndex refers to
    std::array<NeighborType, 27> _neighborTypes{};
    std::array<int, 27> _neighborRanks{};
    std::reference_wrapper<const readdy::model::Context> _context;

public:
    MPIDomain(int rank, int worldSize, std::array<scalar, 3> minDomainWidths, const readdy::model::Context &ctx)
            : _rank(rank), _worldSize(worldSize), _context(std::cref(ctx)), _haloThickness(ctx.calculateMaxCutoff()),
              _nUsedRanks(0) {
        const auto &boxSize = _context.get().boxSize();
        const auto &periodic = _context.get().periodicBoundaryConditions();

        if (_haloThickness <= 0.) {
            // this scenario is anyway pointless because there are no interactions,
            // but we cannot have minDomainWidths = 0
            auto minBoxLength = std::min_element(boxSize.begin(), boxSize.end());
            _haloThickness = *minBoxLength;
        }

        for (int i = 0; i < 3; ++i) {
            if (minDomainWidths[i] <= 0.) {
                minDomainWidths[i] = 2. * _haloThickness;
                readdy::log::info("Setting minDomainWidths[{}] to 2 * haloThickness", i);
            }
        }

        readdy::log::info("minDomainWidths are now {} {} {}",
                minDomainWidths[0], minDomainWidths[1], minDomainWidths[2]);

        if (worldSize >= 2) {
            int coord = 0;

            for (std::size_t i = 0; i < 3; ++i) {
                _nDomainsPerAxis[i] = static_cast<unsigned int>(std::max(1., std::floor(boxSize[i] / minDomainWidths[i])));
            }
            _nUsedRanks = _nDomainsPerAxis[0] * _nDomainsPerAxis[1] * _nDomainsPerAxis[2] + 1;

            while (_nUsedRanks > worldSize) {
                // try to increase size of domains (round robbing), to decrease usedRanks
                minDomainWidths[coord] = 1.5 * minDomainWidths[coord];
                coord = (coord + 1) % 3;
                for (std::size_t i = 0; i < 3; ++i) {
                    _nDomainsPerAxis[i] = static_cast<unsigned int>(std::max(1., std::floor(boxSize[i] / minDomainWidths[i])));
                }
                _nUsedRanks = _nDomainsPerAxis[0] * _nDomainsPerAxis[1] * _nDomainsPerAxis[2] + 1;
            }
            if (not isValidDecomposition(_nDomainsPerAxis)) {
                throw std::runtime_error("Could not determine a valid domain decomposition");
            }

        } else {
            throw std::logic_error("worldSize must be at least 2 (one master, one worker)");
        }

        _nWorkerRanks = _nUsedRanks-1;
        _workerRanks.resize(_nWorkerRanks); // master rank 0 is also used, but not a worker, thus subtract it here
        std::iota(_workerRanks.begin(), _workerRanks.end(), 1);
        assert(_workerRanks[0] == 1);
        assert(_workerRanks.back() == _nUsedRanks-1);
        _nIdleRanks = worldSize - _nUsedRanks;

        _domainIndex = readdy::util::Index3D(_nDomainsPerAxis[0], _nDomainsPerAxis[1], _nDomainsPerAxis[2]);

        // the rest is only for workers
        if (rank != 0 and rank < _nUsedRanks) {
            // find out which this ranks' ijk coordinates are, consider -1 because of master rank 0
            _myIdx = _domainIndex.inverse(rank - 1);
            for (std::size_t i = 0; i < 3; ++i) {
                _extent[i] = boxSize[i] / static_cast<scalar>(_nDomainsPerAxis[i]);
                _origin[i] = -0.5 * boxSize[i] + _myIdx[i] * _extent[i];
                _originWithHalo[i] = _origin[i] - _haloThickness;
                _extentWithHalo[i] = _extent[i] + 2 * _haloThickness;
            }

            // set up neighbors, i.e. the adjacency between domains
            for (int di = -1; di < 2; ++di) {
                for (int dj = -1; dj < 2; ++dj) {
                    for (int dk = -1; dk < 2; ++dk) {
                        int i = _myIdx[0] + di;
                        int j = _myIdx[1] + dj;
                        int k = _myIdx[2] + dk;

                        i = wrapDomainIdx(i, 0);
                        j = wrapDomainIdx(j, 1);
                        k = wrapDomainIdx(k, 2);

                        // determine if neighbor is to be considered as such
                        int otherRank;
                        NeighborType neighborType;
                        if (i == -1 or j == -1 or k == -1) {
                            // other domain is not a neighbor
                            otherRank = -1;
                            neighborType = NeighborType::nan;
                        } else {
                            otherRank = _domainIndex(i, j, k) + 1; // +1 considers master rank
                            if (otherRank == rank) {
                                neighborType = NeighborType::self;
                            } else {
                                neighborType = NeighborType::regular;
                            }
                        }

                        auto dijk = neighborIndex(di + 1, dj + 1, dk + 1);
                        _neighborRanks.at(dijk) = otherRank;
                        _neighborTypes.at(dijk) = neighborType;
                    }
                }
            }

            assert(_neighborRanks.at(neighborIndex(1, 1, 1)) == rank);

        } else if (rank == 0) {
            // master rank 0 must at least know how big domains are
            for (std::size_t i = 0; i < 3; ++i) {
                _extent[i] = boxSize[i] / static_cast<scalar>(_nDomainsPerAxis[i]);
            }
        } else {
            // allocated but unneeded workers
            _isNotIdle = false;
        }
    }

    int rankOfPosition(const Vec3 &pos) const {
        const auto ijk = ijkOfPosition(pos);
        return _domainIndex(ijk[0], ijk[1], ijk[2]) + 1; // + 1 because master rank = 0
    };

    bool isInDomainCore(const Vec3 &pos) const {
        validateRankNotMaster();
        return (_origin.x <= pos.x and pos.x < _origin.x + _extent.x and
                _origin.y <= pos.y and pos.y < _origin.y + _extent.y and
                _origin.z <= pos.z and pos.z < _origin.z + _extent.z);
    }

    // @todo consider hyperplanes
    bool isInDomainCoreOrHalo(const Vec3 &pos) const {
        validateRankNotMaster();
        // additionally need to consider wrapped position if it is at the edge of the box
        // (i.e. in the outer shell with haloThickness of the box)
        // this wrapping however must only be attempted in axes where the domainIdx[axis] of
        // domain and position is non-identical, otherwise the position might be wrapped somewhere else

        // first attempt normal position, which in most cases is sufficient
        if (_originWithHalo.x <= pos.x and pos.x < _originWithHalo.x + _extentWithHalo.x and
             _originWithHalo.y <= pos.y and pos.y < _originWithHalo.y + _extentWithHalo.y and
             _originWithHalo.z <= pos.z and pos.z < _originWithHalo.z + _extentWithHalo.z) {
            return true;
        }
        auto wrappedPos = wrapIntoThisHalo(pos);
        return (_originWithHalo.x <= wrappedPos.x and wrappedPos.x < _originWithHalo.x + _extentWithHalo.x and
                _originWithHalo.y <= wrappedPos.y and wrappedPos.y < _originWithHalo.y + _extentWithHalo.y and
                _originWithHalo.z <= wrappedPos.z and wrappedPos.z < _originWithHalo.z + _extentWithHalo.z);
    }

    Vec3 wrapIntoThisHalo(const Vec3 &pos) const {
        validateRankNotMaster();
        Vec3 wrappedPos(pos);
        const auto &box = _context.get().boxSize();
        const auto &periodic = _context.get().periodicBoundaryConditions();
        const auto ijk = ijkOfPosition(pos);
        for (int d = 0; d < 3; ++d) {
            if (ijk[d] != _myIdx[d]) { // only attempt wrap in axes, which are not identical
                if (periodic[d]) {
                    if (pos[d] < - 0.5 * box[d] + _haloThickness) { // account for halo region in target direction
                        wrappedPos[d] += box[d];
                    } else if (pos[d] > 0.5 * box[d] - _haloThickness) {
                        wrappedPos[d] -= box[d];
                    }
                }
            }
        }
        return wrappedPos;
    }

    [[nodiscard]] bool isInDomainHalo(const Vec3 &pos) const {
        return isInDomainCoreOrHalo(pos) and not isInDomainCore(pos);
    }

    [[nodiscard]] const Vec3 &origin() const {
        validateRankNotMaster();
        return _origin;
    }

    [[nodiscard]] const Vec3 &extent() const {
        // master rank has this information
        return _extent;
    }

    [[nodiscard]] const Vec3 &originWithHalo() const {
        validateRankNotMaster();
        return _originWithHalo;
    }

    [[nodiscard]] const Vec3 &extentWithHalo() const {
        validateRankNotMaster();
        return _extentWithHalo;
    }

    [[nodiscard]] const readdy::util::Index3D &domainIndex() const {
        return _domainIndex;
    }

    [[nodiscard]] const std::array<std::size_t, 3> &nDomainsPerAxis() const {
        return _nDomainsPerAxis;
    }

    [[nodiscard]] const std::size_t nDomains() const {
        return std::accumulate(_nDomainsPerAxis.begin(), _nDomainsPerAxis.end(), 1, std::multiplies<>());
    }

    [[nodiscard]] int nUsedRanks() const {
        return _nUsedRanks;
    }

    [[nodiscard]] int nWorkerRanks() const {
        return _nWorkerRanks;
    }

    [[nodiscard]] int nIdleRanks() const {
        return _nIdleRanks;
    }

    [[nodiscard]] bool isMasterRank() const {
        return (_rank == 0);
    }

    [[nodiscard]] bool isWorkerRank() const {
        return (_isNotIdle and not isMasterRank());
    }

    [[nodiscard]] bool isIdleRank() const {
        return (not _isNotIdle);
    }

    [[nodiscard]] const std::vector<int> &workerRanks() const {
        return _workerRanks;
    }

    /** calculate the core region (given by origin and extent) of domain associated with otherRank */
    [[nodiscard]] std::pair<Vec3, Vec3> coreOfDomain(int otherRank) const {
        if (otherRank != 0 and otherRank < _nUsedRanks) {
            // find out which this ranks' ijk coordinates are, consider -1 because of master rank 0
            const auto &boxSize = _context.get().boxSize();
            auto ijkOfOtherRank = _domainIndex.inverse(otherRank - 1);
            Vec3 origin, extent;
            for (std::size_t i = 0; i < 3; ++i) {
                extent[i] = boxSize[i] / static_cast<scalar>(_nDomainsPerAxis[i]);
                origin[i] = -0.5 * boxSize[i] + ijkOfOtherRank[i] * extent[i];
                //originWithHalo[i] = origin[i] - haloThickness;
                //extentWithHalo[i] = extent[i] + 2 * haloThickness;
            }
            return std::make_pair(origin, extent);
        } else {
            throw std::runtime_error(fmt::format("Can only determine core region of a worker rank"));
        }
    }

    /**
     * @return This domains' (ijk) array
     */
    [[nodiscard]] const std::array<std::size_t, 3> &myIdx() const {
        validateRankNotMaster();
        return _myIdx;
    }

    [[nodiscard]] const std::array<NeighborType, 27> &neighborTypes() const {
        validateRankNotMaster();
        return _neighborTypes;
    }

    [[nodiscard]] const std::array<int, 27> &neighborRanks() const {
        validateRankNotMaster();
        return _neighborRanks;
    }

    [[nodiscard]] std::array<std::size_t, 3> ijkOfPosition(const Vec3 &pos) const {
        const auto &boxSize = _context.get().boxSize();
        if (!(-.5 * boxSize[0] <= pos.x && .5 * boxSize[0] > pos.x
              && -.5 * boxSize[1] <= pos.y && .5 * boxSize[1] > pos.y
              && -.5 * boxSize[2] <= pos.z && .5 * boxSize[2] > pos.z)) {
            throw std::logic_error(fmt::format("ijkOfPosition: position {} was out of bounds.", pos));
        }
        return {
                        static_cast<std::size_t>(std::floor((pos.x + .5 * boxSize[0]) / _extent.x)),
                        static_cast<std::size_t>(std::floor((pos.y + .5 * boxSize[1]) / _extent.y)),
                        static_cast<std::size_t>(std::floor((pos.z + .5 * boxSize[2]) / _extent.z))
                };
    }

private:
    void validateRankNotMaster() const {
        if (_rank == 0) {
            throw std::logic_error("Master rank 0 cannot know which domain you're referring to.");
        }
    }

    [[nodiscard]] int wrapDomainIdx(int posIdx, std::uint8_t axis) const {
        auto &pbc = _context.get().periodicBoundaryConditions();
        auto nDomainsAxis = static_cast<int>(_domainIndex[axis]);
        if (pbc[axis]) {
            return (posIdx % nDomainsAxis + nDomainsAxis) % nDomainsAxis;
        } else if (posIdx < 0 or posIdx >= nDomainsAxis) {
            return -1;
        } else {
            return posIdx;
        }
    }

    [[nodiscard]] bool isValidDecomposition(const std::array<std::size_t, 3> nDomains) const {
        const auto cutoff = _context.get().calculateMaxCutoff();
        const auto periodic = _context.get().periodicBoundaryConditions();
        const auto boxSize = _context.get().boxSize();
        std::array<scalar, 3> domainWidths{};
        for (std::size_t i = 0; i < 3; ++i) {
            domainWidths[i] = boxSize[i] / static_cast<scalar>(nDomains[i]);
            if (domainWidths[i] < 2. * cutoff) { // for halo regions to make sense and not overlap
                if (nDomains[i] == 1 and not periodic[i]) {
                    /* smaller than cutoff is ok, when there are no neighbors to be considered */
                } else {
                    return false;
                }
            }
        }
        return true;
    }
};

}
