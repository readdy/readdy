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
 * @file MPIDomain.h
 * @brief « brief description »
 * @author chrisfroe
 * @date 17.06.19
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/model/Context.h>

namespace readdy::kernel::mpi::model {

class MPIDomain {
public:

    MPIDomain(int rank, int worldSize, std::array<scalar, 3> minDomainWidths, const readdy::model::Context &ctx)
            : rank(rank), worldSize(worldSize), _context(std::cref(ctx)), haloThickness(ctx.calculateMaxCutoff()) {
        const auto &boxSize = _context.get().boxSize();
        const auto &periodic = _context.get().periodicBoundaryConditions();

        // Find out nDomains per axis from user given widths
        {
            std::array<scalar, 3> domainWidths{};
            for (std::size_t i = 0; i < 3; ++i) {
                nDomains[i] = static_cast<unsigned int>(std::max(1., std::floor(boxSize[i] / minDomainWidths[i])));
                domainWidths[i] = boxSize[i] / static_cast<scalar>(nDomains[i]);
            }

            if (rank == 0) {
                readdy::log::info("MPI spatial domain decomposition:");
                readdy::log::info("The user given minimal domain widths were ({}, {}, {})",
                                  minDomainWidths[0],
                                  minDomainWidths[1],
                                  minDomainWidths[2]);
                readdy::log::info("there will be {} * {} * {} = {} number of domains", nDomains[0], nDomains[1],
                                  nDomains[2], nDomains[0] * nDomains[1] * nDomains[2]);
                readdy::log::info("with actual widths dx {} dy {} dz {}", domainWidths[0], domainWidths[1],
                                  domainWidths[2]);
            }

            if (not isValidDecomposition(nDomains)) {
                throw std::logic_error("Spatial decomposition is not valid.");
            }

            const auto numberDomains = nDomains[0] * nDomains[1] * nDomains[2];
            if (numberDomains + 1 != worldSize) {// add one for master rank 0
                throw std::logic_error(
                        fmt::format("There are {} + 1 worker positions to be filled, but there are {} workers",
                                    numberDomains + 1, worldSize));
            }
        }

        domainIndex = util::Index3D(nDomains[0], nDomains[1], nDomains[2]);

        // the rest is only for workers
        if (rank != 0) {
            // find out which this ranks' ijk coordinates are, consider -1 because of master rank 0
            _myIdx = domainIndex.inverse(rank - 1);
            for (std::size_t i = 0; i < 3; ++i) {
                _extent[i] = boxSize[i] / static_cast<scalar>(nDomains[i]);
                _origin[i] = -0.5 * boxSize[i] + _myIdx[i] * _extent[i];
                _originWithHalo[i] = _origin[i] - haloThickness;
                _extentWithHalo[i] = _extent[i] + 2 * haloThickness;
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
                            otherRank = domainIndex(i, j, k) + 1; // +1 considers master rank
                            if (otherRank == rank) {
                                neighborType = NeighborType::self;
                            } else {
                                neighborType = NeighborType::regular;
                            }
                        }

                        auto dijk = _neighborIndex(di + 1, dj + 1, dk + 1);
                        _neighborRanks.at(dijk) = otherRank;
                        _neighborTypes.at(dijk) = neighborType;
                    }
                }
            }

            assert(_neighborRanks.at(_neighborIndex(1, 1, 1)) == rank);

        } else {
            // master rank 0 must at least know how big domains are
            for (std::size_t i = 0; i < 3; ++i) {
                _extent[i] = boxSize[i] / static_cast<scalar>(nDomains[i]);
            }
        }
    }

    int rankOfPosition(const Vec3 &pos) const {
        const auto &boxSize = _context.get().boxSize();
        if (!(-.5 * boxSize[0] <= pos.x && .5 * boxSize[0] > pos.x
              && -.5 * boxSize[1] <= pos.y && .5 * boxSize[1] > pos.y
              && -.5 * boxSize[2] <= pos.z && .5 * boxSize[2] > pos.z)) {
            throw std::logic_error("rankOfPosition: position was out of bounds.");
        }
        const auto i = static_cast<std::size_t>(std::floor((pos.x + .5 * boxSize[0]) / _extent.x));
        const auto j = static_cast<std::size_t>(std::floor((pos.y + .5 * boxSize[1]) / _extent.y));
        const auto k = static_cast<std::size_t>(std::floor((pos.z + .5 * boxSize[2]) / _extent.z));
        return domainIndex(i, j, k);
    };

    bool isInDomainCore(const Vec3 &pos) const {
        if (rank != 0) {
            return (_origin.x <= pos.x and pos.x < _origin.x + _extent.x and
                    _origin.y <= pos.y and pos.y < _origin.y + _extent.y and
                    _origin.z <= pos.z and pos.z < _origin.z + _extent.z);
        } else {
            throw std::logic_error("Master rank 0 cannot know which domain you're referring to.");
        }
    }

    bool isInDomainCoreOrHalo(const Vec3 &pos) const {
        if (rank != 0) {
            return (_originWithHalo.x <= pos.x and pos.x < _originWithHalo.x + _extentWithHalo.x and
                    _originWithHalo.y <= pos.y and pos.y < _originWithHalo.y + _extentWithHalo.y and
                    _originWithHalo.z <= pos.z and pos.z < _originWithHalo.z + _extentWithHalo.z);
        } else {
            throw std::logic_error("Master rank 0 cannot know which domain you're referring to.");
        }
    }

    bool isInDomainHalo(const Vec3 &pos) {
        return isInDomainCoreOrHalo(pos) and not isInDomainCore(pos);
    }

    const int rank;
    const int worldSize;
    const scalar haloThickness;
    std::array<std::size_t, 3> nDomains{};
    util::Index3D domainIndex; // rank of (ijk) is domainIndex(i,j,k)+1

private:
    int wrapDomainIdx(int posIdx, std::uint8_t axis) const {
        auto &pbc = _context.get().periodicBoundaryConditions();
        auto nDomainsAxis = static_cast<int>(domainIndex[axis]);
        if (pbc[axis]) {
            return (posIdx % nDomainsAxis + nDomainsAxis) % nDomainsAxis;
        } else if (posIdx < 0 or posIdx >= nDomainsAxis) {
            return -1;
        } else {
            return posIdx;
        }
    }

    bool isValidDecomposition(const std::array<std::size_t, 3> dims) const {
        const auto cutoff = _context.get().calculateMaxCutoff();
        const auto periodic = _context.get().periodicBoundaryConditions();
        const auto boxSize = _context.get().boxSize();
        std::array<scalar, 3> boxWidths{};
        for (std::size_t i = 0; i < 3; ++i) {
            boxWidths[i] = boxSize[i] / static_cast<scalar>(dims[i]);
            if (boxWidths[i] < 2. * cutoff) { // for halo regions to make sense and not overlap
                if (dims[i] == 1 and not periodic[0]) {
                    /* smaller than cutoff is ok, when there are no neighbors to be considered */
                } else {
                    return false;
                }
            }
        }
        return true;
    }

    /** The following will only be defined for rank != 0 */
    // origin and extent define the core region of the domain
    Vec3 _origin; // lower-left corner
    Vec3 _extent;
    // these define the extended region core+halo
    Vec3 _originWithHalo;
    Vec3 _extentWithHalo;

    std::array<std::size_t, 3> _myIdx{}; // (ijk) indices of this domain
    const std::size_t _three{3};

    // map from ([0-2], [0-2], [0-2]) to the index of the 27 neighbors (including self)
    const util::Index3D _neighborIndex = util::Index3D(_three, _three, _three);
    enum NeighborType {
        self, // is a valid neighbor, but due to periodicity it is this domain
        nan, // not a neighbor (i.e. end of simulation box and no periodic bound.)
        regular // neighbor is another domain ()
    };
    std::array<NeighborType, 27> _neighborTypes{};
    std::array<int, 27> _neighborRanks{};
    std::reference_wrapper<const readdy::model::Context> _context;
};

}