/**
 * << detailed description >>
 *
 * @file NeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.07.16
 */

#ifndef READDY_MAIN_NEIGHBORLIST_H
#define READDY_MAIN_NEIGHBORLIST_H

#include <memory>
#include <readdy/common/make_unique.h>
#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/kernel/cpu/model/ParticleData.h>
#include <readdy/model/KernelContext.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace model {
                template<typename container=std::vector<ParticleIndexPair>>
                struct NeighborListBase {
                    typedef typename container::iterator iter_type;
                    typedef typename container::const_iterator const_iter_type;

                    virtual iter_type begin() { return pairs->begin(); };
                    virtual const_iter_type begin() const { return cbegin(); };
                    virtual const_iter_type cbegin() const { return pairs->cbegin(); };
                    virtual iter_type end() { return pairs->end(); };
                    virtual const_iter_type end() const { return cend(); };
                    virtual const_iter_type cend() const { return pairs->cend(); };

                    virtual void create(const ParticleData &data) = 0;
                protected:
                    std::unique_ptr<container> pairs = std::make_unique<container>();
                };

                class NeighborList : public NeighborListBase<> {
                public:
                    NeighborList(readdy::model::KernelContext *const context) : ctx(context) {
                    }

                    virtual void create(const ParticleData &data) override {
                        setupBoxes();
                        fillBoxes(data);
                    }

                    virtual void setupBoxes() {
                        const auto simBoxSize = ctx->getBoxSize();
                        if (boxes.empty()) {
                            double maxCutoff = 0;
                            for (auto &&e : ctx->getAllOrder2RegisteredPotentialTypes()) {
                                for (auto &&p : ctx->getOrder2Potentials(std::get<0>(e), std::get<1>(e))) {
                                    maxCutoff = maxCutoff < p->getCutoffRadius() ? p->getCutoffRadius() : maxCutoff;
                                }
                            }
                            for (auto &&e : ctx->getAllOrder2Reactions()) {
                                maxCutoff = maxCutoff < e->getEductDistance() ? e->getEductDistance() : maxCutoff;
                            }
                            NeighborList::maxCutoff = maxCutoff;
                            if (maxCutoff > 0) {

                                for (unsigned int i = 0; i < 3; ++i) {
                                    nBoxes[i] = (int) floor(simBoxSize[i] / maxCutoff);
                                    if (nBoxes[i] == 0) nBoxes[i] = 1;
                                    boxSize[i] = simBoxSize[i] / nBoxes[i];
                                }
                                for (long i = 0; i < nBoxes[0]; ++i) {
                                    for (long j = 0; j < nBoxes[1]; ++j) {
                                        for (long k = 0; k < nBoxes[2]; ++k) {
                                            boxes.push_back({i, j, k, k + j * nBoxes[2] +
                                                                      i * nBoxes[2] * nBoxes[1]});
                                        }
                                    }
                                }
                                for (long i = 0; i < nBoxes[0]; ++i) {
                                    for (long j = 0; j < nBoxes[1]; ++j) {
                                        for (long k = 0; k < nBoxes[2]; ++k) {
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
                                    }
                                }
                            }
                        }
                    }
                    virtual void fillBoxes(const ParticleData &data) {
                        const auto simBoxSize = ctx->getBoxSize();
                        if (maxCutoff > 0) {

                            for (auto &&box : boxes) {
                                box.particleIndices.clear();
                            }

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
                                        pairs->push_back({pI, box.particleIndices[j]});
                                    }

                                    for (auto &&neighboringBox : box.neighboringBoxes) {
                                        for (const auto &pJ : neighboringBox->particleIndices) {
                                            pairs->push_back({pI, pJ});
                                        }
                                    }
                                }
                            }
                        }
                    }


                protected:
                    struct Box {
                        std::vector<Box *> neighboringBoxes{};
                        std::vector<unsigned long> particleIndices{};
                        long i, j, k;
                        long id = 0;

                        Box(long i, long j, long k, long id) : i(i), j(j), k(k), id(id) {
                        }
                        void addNeighbor(Box *box) {
                            if (box && box->id != id) neighboringBoxes.push_back(box);
                        }
                    };
                    std::vector<Box> boxes{};
                    std::array<int, 3> nBoxes{{0, 0, 0}};
                    readdy::model::Vec3 boxSize{0, 0, 0};
                    double maxCutoff = 0;

                    long positive_modulo(long i, long n) const {
                        return (i % n + n) % n;
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
                    readdy::model::KernelContext *const ctx;
                };
            }
        }
    }
}
#endif //READDY_MAIN_NEIGHBORLIST_H
