/**
 * << detailed description >>
 *
 * @file SingleCPUNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#ifndef READDY_MAIN_SINGLECPUNEIGHBORLIST_H
#define READDY_MAIN_SINGLECPUNEIGHBORLIST_H

#include "SingleCPUParticleData.h"
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <readdy/model/KernelContext.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace model {
struct ParticleIndexPair {
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
        boost::hash_combine(seed, pip.idx1);
        boost::hash_combine(seed, pip.idx2);
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

struct ParticleIndexPairHasher {
    size_t operator()(const ParticleIndexPair &pip) const {
        return hash_value(pip);
    }
};

template<typename container=std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher>>
struct SingleCPUNeighborListContainer {
    typedef typename container::iterator iter_type;
    typedef typename container::const_iterator const_iter_type;

    virtual iter_type begin() { return pairs->begin(); };

    virtual const_iter_type begin() const { return cbegin(); };

    virtual const_iter_type cbegin() const { return pairs->cbegin(); };

    virtual iter_type end() { return pairs->end(); };

    virtual const_iter_type end() const { return cend(); };

    virtual const_iter_type cend() const { return pairs->cend(); };

    virtual void create(const SingleCPUParticleData &data) = 0;

protected:
    std::unique_ptr<container> pairs = std::make_unique<container>();
};

struct NaiveSingleCPUNeighborList : public SingleCPUNeighborListContainer<> {
    virtual void create(const SingleCPUParticleData &data) override;

    // ctor and dtor
    NaiveSingleCPUNeighborList();

    ~NaiveSingleCPUNeighborList();

    // move
    NaiveSingleCPUNeighborList(NaiveSingleCPUNeighborList &&rhs);

    NaiveSingleCPUNeighborList &operator=(NaiveSingleCPUNeighborList &&rhs);

    // copy
    NaiveSingleCPUNeighborList(const NaiveSingleCPUNeighborList &rhs) = delete;

    NaiveSingleCPUNeighborList &operator=(const NaiveSingleCPUNeighborList &rhs) = delete;

};

struct Box {
    std::vector<Box *> neighbors{};
    std::vector<unsigned long> particleIndices{};
    long i, j, k;
    long id = 0;

    Box(long i, long j, long k, long id) : i(i), j(j), k(k), id(id) {
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
class NotThatNaiveSingleCPUNeighborList : public SingleCPUNeighborListContainer<container> {
    using super = readdy::kernel::singlecpu::model::SingleCPUNeighborListContainer<container>;
    using context_t = readdy::model::KernelContext;
public:
    NotThatNaiveSingleCPUNeighborList(const context_t *const context) : ctx(context) {
    }

    virtual void create(const SingleCPUParticleData &data) override {
        setupBoxes();
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

    virtual void setupBoxes() {
        if (boxes.empty()) {
            const auto simBoxSize = ctx->getBoxSize();
            double maxCutoff = 0;
            for (auto &&p : ctx->getVectorAllOrder2Potentials()) {
                maxCutoff = maxCutoff < p->getCutoffRadius() ? p->getCutoffRadius() : maxCutoff;
            }
            for (auto &&e : ctx->getAllOrder2Reactions()) {
                maxCutoff = maxCutoff < e->getEductDistance() ? e->getEductDistance() : maxCutoff;
            }
            NotThatNaiveSingleCPUNeighborList::maxCutoff = maxCutoff;
            if (maxCutoff > 0) {

                for (unsigned int i = 0; i < 3; ++i) {
                    nBoxes[i] = (int) floor(simBoxSize[i] / maxCutoff);
                    if (nBoxes[i] == 0) nBoxes[i] = 1;
                    boxSize[i] = simBoxSize[i] / nBoxes[i];
                }
                for (long i = 0; i < nBoxes[0]; ++i) {
                    for (long j = 0; j < nBoxes[1]; ++j) {
                        for (long k = 0; k < nBoxes[2]; ++k) {
                            boxes.push_back({i, j, k, k + j * nBoxes[2] + i * nBoxes[2] * nBoxes[1]});
                        }
                    }
                }
                for (unsigned long i = 0; i < nBoxes[0]; ++i) {
                    for (unsigned long j = 0; j < nBoxes[1]; ++j) {
                        for (unsigned long k = 0; k < nBoxes[2]; ++k) {
                            setupNeighboringBoxes(i, j, k);
                        }
                    }
                }
            }
        }
    }

    virtual void fillBoxes(const SingleCPUParticleData &data) {
        const auto simBoxSize = ctx->getBoxSize();
        if (maxCutoff > 0) {

            for (auto &&box : boxes) {
                box.particleIndices.clear();
            }
            super::pairs->clear();

            auto it_pos = data.cbegin_positions();
            unsigned long idx = 0;
            const auto shift = readdy::model::Vec3(.5 * simBoxSize[0], .5 * simBoxSize[1],
                                                   .5 * simBoxSize[2]);
            while (it_pos != data.cend_positions()) {
                const auto pos_shifted = *it_pos + shift;
                const long i = (const long) floor(pos_shifted[0] / boxSize[0]);
                const long j = (const long) floor(pos_shifted[1] / boxSize[1]);
                const long k = (const long) floor(pos_shifted[2] / boxSize[2]);
                auto box = getBox(i, j, k);
                if (box) {
                    box->particleIndices.push_back(idx);
                }
                ++idx;
                ++it_pos;
            }

            for (auto &&box : boxes) {
                for (long i = 0; i < box.particleIndices.size(); ++i) {
                    const auto pI = box.particleIndices[i];
                    for (long j = i + 1; j < box.particleIndices.size(); ++j) {
                        super::pairs->insert((*super::pairs).end(), {pI, box.particleIndices[j]});
                    }

                    for (auto &&neighboringBox : box.neighbors) {
                        for (const auto &pJ : neighboringBox->particleIndices) {
                            super::pairs->insert((*super::pairs).end(), {pI, pJ});
                        }
                    }
                }
            }
        }
    }

    Box *getBox(long i, long j, long k) {
        const auto &periodic = ctx->getPeriodicBoundary();
        if (periodic[0]) i = positive_modulo(i, nBoxes[0]);
        else if (i < 0 || i >= nBoxes[0]) return nullptr;
        if (periodic[1]) j = positive_modulo(j, nBoxes[1]);
        else if (j < 0 || j >= nBoxes[1]) return nullptr;
        if (periodic[2]) k = positive_modulo(k, nBoxes[2]);
        else if (k < 0 || k >= nBoxes[2]) return nullptr;
        return &boxes[k + j * nBoxes[2] + i * nBoxes[2] * nBoxes[1]];
    }

    const std::vector<Box> &getBoxes() const {
        return boxes;
    }

protected:
    std::vector<Box> boxes{};
    std::array<int, 3> nBoxes{{0, 0, 0}};
    readdy::model::Vec3 boxSize{0, 0, 0};
    double maxCutoff = 0;

    long positive_modulo(long i, long n) const {
        return (i % n + n) % n;
    }


    const context_t *const ctx;
};

struct SingleCPUNeighborList : public NotThatNaiveSingleCPUNeighborList<std::vector<ParticleIndexPair>> {
    SingleCPUNeighborList(const readdy::model::KernelContext *const ctx)
            : NotThatNaiveSingleCPUNeighborList(ctx) {}
};

}
}
}
}

#endif //READDY_MAIN_SINGLECPUNEIGHBORLIST_H
