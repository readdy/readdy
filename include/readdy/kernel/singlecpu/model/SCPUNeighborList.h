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
 * << detailed description >>
 *
 * @file SingleCPUNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#pragma once
#include <unordered_set>
#include <readdy/model/Context.h>
#include "SCPUParticleData.h"
#include <readdy/common/numeric.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {
struct READDY_API ParticleIndexPair {
    size_t idx1, idx2;

    ParticleIndexPair(size_t idx1, size_t idx2) {
        if (idx1 < idx2) {
            ParticleIndexPair::idx1 = idx1;
            ParticleIndexPair::idx2 = idx2;
        } else if (idx1 > idx2) {
            ParticleIndexPair::idx1 = idx2;
            ParticleIndexPair::idx2 = idx1;
        } else {
            throw std::runtime_error("pair must not have equal indices");
        }
    }

    friend size_t hash_value(const ParticleIndexPair &pip) {
        size_t seed = 0;
        readdy::util::hash::combine(seed, pip.idx1);
        readdy::util::hash::combine(seed, pip.idx2);
        return seed;
    }

    friend bool operator==(const ParticleIndexPair &pip1, const ParticleIndexPair &pip2) {
        return pip1.idx1 == pip2.idx1 && pip1.idx2 == pip2.idx2;
    }

    friend std::ostream &operator<<(std::ostream &os, const ParticleIndexPair &pip) {
        os << "ParticleIndexPair(" << pip.idx1 << ", " << pip.idx2 << ")";
        return os;
    }
};

struct READDY_API ParticleIndexPairHasher {
    size_t operator()(const ParticleIndexPair &pip) const {
        return hash_value(pip);
    }
};

template<typename container=std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher>>
struct READDY_API SCPUNeighborListContainer {
    typedef typename container::iterator iter_type;
    typedef typename container::const_iterator const_iter_type;

    virtual iter_type begin() { return pairs->begin(); };

    virtual const_iter_type begin() const { return cbegin(); };

    virtual const_iter_type cbegin() const { return pairs->cbegin(); };

    virtual iter_type end() { return pairs->end(); };

    virtual const_iter_type end() const { return cend(); };

    virtual const_iter_type cend() const { return pairs->cend(); };

    virtual void create(const SCPUParticleData &data, scalar skin) = 0;

    scalar skin() const {
        return _skin;
    }

protected:
    std::unique_ptr<container> pairs = std::make_unique<container>();

    scalar _skin;
};

struct READDY_API SCPUNaiveNeighborList : public SCPUNeighborListContainer<> {
    void create(const SCPUParticleData &data, scalar skin) override;

    // ctor and dtor
    SCPUNaiveNeighborList();

    ~SCPUNaiveNeighborList();

    // move
    SCPUNaiveNeighborList(SCPUNaiveNeighborList &&rhs) noexcept;

    SCPUNaiveNeighborList &operator=(SCPUNaiveNeighborList &&rhs) noexcept;

    // copy
    SCPUNaiveNeighborList(const SCPUNaiveNeighborList &rhs) = delete;

    SCPUNaiveNeighborList &operator=(const SCPUNaiveNeighborList &rhs) = delete;

};

struct READDY_API Box {
    std::vector<Box *> neighbors;
    std::vector<unsigned long> particleIndices;
    long i, j, k;
    long id = 0;

    Box(long i, long j, long k, long id) : i(i), j(j), k(k), id(id), particleIndices(), neighbors() {
    }

    void addNeighbor(Box *box) {
        if (box && box->id != id
            && std::find(neighbors.begin(), neighbors.end(), box) == neighbors.end()) {
            neighbors.push_back(box);
        }
    }

    friend bool operator==(const Box &lhs, const Box &rhs) {
        return lhs.id == rhs.id;
    }

    friend bool operator!=(const Box &lhs, const Box &rhs) {
        return !(lhs == rhs);
    }
};

template<typename container=std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher>>
class READDY_API SCPUNotThatNaiveNeighborList : public SCPUNeighborListContainer<container> {
    using super = readdy::kernel::scpu::model::SCPUNeighborListContainer<container>;
    using context = readdy::model::Context;
public:
    explicit SCPUNotThatNaiveNeighborList(const context *const context) : ctx(context) {}

    void create(const SCPUParticleData &data, scalar skin) override {
        super::_skin = skin;
        setupBoxes(skin);
        fillBoxes(data);
    }

    void clear() {
        boxes.clear();
    }

    virtual void setupNeighboringBoxes(unsigned long i, unsigned long j, unsigned long k) {
        auto me = getBox(i, j, k);
        me->addNeighbor(getBox(i + 0, j + 0, k + 1));
        me->addNeighbor(getBox(i + 0, j + 1, k - 1));
        me->addNeighbor(getBox(i + 0, j + 1, k + 0));
        me->addNeighbor(getBox(i + 0, j + 1, k + 1));
        me->addNeighbor(getBox(i + 1, j - 1, k - 1));
        me->addNeighbor(getBox(i + 1, j - 1, k + 0));
        me->addNeighbor(getBox(i + 1, j - 1, k + 1));
        me->addNeighbor(getBox(i + 1, j + 0, k - 1));
        me->addNeighbor(getBox(i + 1, j + 0, k + 0));
        me->addNeighbor(getBox(i + 1, j + 0, k + 1));
        me->addNeighbor(getBox(i + 1, j + 1, k - 1));
        me->addNeighbor(getBox(i + 1, j + 1, k + 0));
        me->addNeighbor(getBox(i + 1, j + 1, k + 1));
    }

    virtual void setupBoxes(scalar skin) {
        if (boxes.empty()) {
            const auto &simBoxSize = ctx->boxSize();
            auto maxCutoff = ctx->calculateMaxCutoff();
            if (maxCutoff > 0) {
                maxCutoff += skin;
                SCPUNotThatNaiveNeighborList::maxCutoff = maxCutoff;

                for (std::uint8_t i = 0; i < 3; ++i) {
                    nBoxes[i] = std::max(1, static_cast<int>(std::floor(simBoxSize[i] / maxCutoff)));
                    boxSize[i] = simBoxSize[i] / nBoxes[i];
                }
                for (long i = 0; i < nBoxes[0]; ++i) {
                    for (long j = 0; j < nBoxes[1]; ++j) {
                        for (long k = 0; k < nBoxes[2]; ++k) {
                            boxes.push_back({i, j, k, k + j * nBoxes[2] + i * nBoxes[2] * nBoxes[1]});
                        }
                    }
                }
                for (auto i = 0; i < nBoxes[0]; ++i) {
                    for (auto j = 0; j < nBoxes[1]; ++j) {
                        for (auto k = 0; k < nBoxes[2]; ++k) {
                            setupNeighboringBoxes(static_cast<unsigned long>(i),
                                                  static_cast<unsigned long>(j),
                                                  static_cast<unsigned long>(k));
                        }
                    }
                }
            }
        }
    }

    virtual void fillBoxes(const SCPUParticleData &data) {
        const auto &simBoxSize = ctx->boxSize();
        if (maxCutoff > 0) {

            std::for_each(boxes.begin(), boxes.end(), [](Box& box) {box.particleIndices.clear();});
            super::pairs->clear();

            unsigned long idx = 0;
            const auto shift = Vec3(c_::half * simBoxSize[0], c_::half * simBoxSize[1], c_::half * simBoxSize[2]);
            for(const auto& entry : data) {
                if(!entry.is_deactivated()) {
                    const auto pos_shifted = entry.position() + shift;
                    const long i = (const long) floor(pos_shifted[0] / boxSize[0]);
                    const long j = (const long) floor(pos_shifted[1] / boxSize[1]);
                    const long k = (const long) floor(pos_shifted[2] / boxSize[2]);
                    auto box = getBox(i, j, k);
                    if (box) {
                        box->particleIndices.push_back(idx);
                    }
                }
                ++idx;
            }

            for (auto &&box : boxes) {
                for (long i = 0; i < box.particleIndices.size(); ++i) {
                    const auto pI = box.particleIndices[i];
                    for (long j = i + 1; j < box.particleIndices.size(); ++j) {
                        super::pairs->push_back({pI, box.particleIndices[j]});
                    }

                    for (auto &&neighboringBox : box.neighbors) {
                        for (const auto &pJ : neighboringBox->particleIndices) {
                            super::pairs->push_back({pI, pJ});
                        }
                    }
                }
            }
        }
    }

    Box *getBox(long i, long j, long k) {
        const auto &periodic = ctx->periodicBoundaryConditions();
        if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nBoxes[0]);
        else if (i < 0 || i >= nBoxes[0]) return nullptr;
        if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nBoxes[1]);
        else if (j < 0 || j >= nBoxes[1]) return nullptr;
        if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nBoxes[2]);
        else if (k < 0 || k >= nBoxes[2]) return nullptr;
        return &boxes[k + j * nBoxes[2] + i * nBoxes[2] * nBoxes[1]];
    }

    const std::vector<Box> &getBoxes() const {
        return boxes;
    }

protected:
    std::vector<Box> boxes{};
    std::array<int, 3> nBoxes{{0, 0, 0}};
    Vec3 boxSize{0, 0, 0};
    scalar maxCutoff {0};

    const context *const ctx;
};

struct READDY_API SCPUNeighborList : public SCPUNotThatNaiveNeighborList<std::vector<ParticleIndexPair>> {
    explicit SCPUNeighborList(const readdy::model::Context *const ctx)
            : SCPUNotThatNaiveNeighborList(ctx) {}
};

}
}
}
}
