/**
 * << detailed description >>
 *
 * @file SingleCPUNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>
#include <iostream>
#include <readdy/common/make_unique.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace model {

                struct NaiveSingleCPUNeighborList::Impl {
                    std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher> pairs{};
                };

                void NaiveSingleCPUNeighborList::create(const SingleCPUParticleData &data) {
                    pimpl->pairs.clear();
                    for (size_t i = 0; i < data.size(); ++i) {
                        for (size_t j = i + 1; j < data.size(); ++j) {
                            pimpl->pairs.emplace(i, j);
                        }
                    }
                }

                NaiveSingleCPUNeighborList &NaiveSingleCPUNeighborList::operator=(
                        NaiveSingleCPUNeighborList &&rhs) = default;


                NaiveSingleCPUNeighborList::NaiveSingleCPUNeighborList() : pimpl(std::make_unique<Impl>()) { }

                iter_type NaiveSingleCPUNeighborList::begin() {
                    return pimpl->pairs.begin();
                }

                const_iter_type NaiveSingleCPUNeighborList::begin() const {
                    return cbegin();
                }

                const_iter_type NaiveSingleCPUNeighborList::cbegin() const {
                    return pimpl->pairs.cbegin();
                }

                iter_type NaiveSingleCPUNeighborList::end() {
                    return pimpl->pairs.end();
                }

                const_iter_type NaiveSingleCPUNeighborList::end() const {
                    return cend();
                }

                const_iter_type NaiveSingleCPUNeighborList::cend() const {
                    return pimpl->pairs.cend();
                }

                NaiveSingleCPUNeighborList::~NaiveSingleCPUNeighborList() = default;

                NaiveSingleCPUNeighborList::NaiveSingleCPUNeighborList(NaiveSingleCPUNeighborList &&rhs) = default;

                struct NotThatNaiveSingleCPUNeighborList::Box {

                    std::vector<Box *> neighboringBoxes{};
                    std::vector<long> particleIndices{};
                    long i, j, k;
                    long id = 0;

                    Box(long i, long j, long k, long id) : i(i), j(j), k(k), id(id) {
                    };

                    void addNeighbor(Box *box) {
                        if (box && box->id != id) neighboringBoxes.push_back(box);
                    }
                };

                struct NotThatNaiveSingleCPUNeighborList::Impl {
                    const readdy::model::KernelContext *ctx;
                    std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher> pairs{};
                    std::vector<Box> boxes{};
                    std::array<int, 3> nBoxes{{0, 0, 0}};
                    readdy::model::Vec3 boxSize{0, 0, 0};
                    double maxCutoff = 0;


                    inline long positive_modulo(long i, long n) const {
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
                };

                NotThatNaiveSingleCPUNeighborList::NotThatNaiveSingleCPUNeighborList(
                        const readdy::model::KernelContext *const ctx) : pimpl(std::make_unique<Impl>()) {
                    pimpl->ctx = ctx;
                }

                void NotThatNaiveSingleCPUNeighborList::create(const SingleCPUParticleData &data) {
                    const auto simBoxSize = pimpl->ctx->getBoxSize();
                    if (pimpl->boxes.empty()) {
                        double maxCutoff = 0;
                        for (auto &&e : pimpl->ctx->getAllOrder2RegisteredPotentialTypes()) {
                            for (auto &&p : pimpl->ctx->getOrder2Potentials(std::get<0>(e), std::get<1>(e))) {
                                maxCutoff = maxCutoff < p->getCutoffRadius() ? p->getCutoffRadius() : maxCutoff;
                            }
                        }
                        for (auto &&e : pimpl->ctx->getAllOrder2Reactions()) {
                            maxCutoff = maxCutoff < e->getEductDistance() ? e->getEductDistance() : maxCutoff;
                        }
                        pimpl->maxCutoff = maxCutoff;
                        if (maxCutoff > 0) {

                            for (unsigned int i = 0; i < 3; ++i) {
                                pimpl->nBoxes[i] = (int) floor(simBoxSize[i] / maxCutoff);
                                if (pimpl->nBoxes[i] == 0) pimpl->nBoxes[i] = 1;
                                pimpl->boxSize[i] = simBoxSize[i] / pimpl->nBoxes[i];
                            }
                            for (long i = 0; i < pimpl->nBoxes[0]; ++i) {
                                for (long j = 0; j < pimpl->nBoxes[1]; ++j) {
                                    for (long k = 0; k < pimpl->nBoxes[2]; ++k) {
                                        pimpl->boxes.push_back({i, j, k, k + j * pimpl->nBoxes[2] +
                                                                         i * pimpl->nBoxes[2] * pimpl->nBoxes[1]});
                                    }
                                }
                            }
                            for (long i = 0; i < pimpl->nBoxes[0]; ++i) {
                                for (long j = 0; j < pimpl->nBoxes[1]; ++j) {
                                    for (long k = 0; k < pimpl->nBoxes[2]; ++k) {
                                        auto me = pimpl->getBox(i, j, k);
                                        me->addNeighbor(pimpl->getBox(i + 0, j + 0, k + 1));
                                        me->addNeighbor(pimpl->getBox(i + 0, j + 1, k - 1));
                                        me->addNeighbor(pimpl->getBox(i + 0, j + 1, k + 0));
                                        me->addNeighbor(pimpl->getBox(i + 0, j + 1, k + 1));
                                        me->addNeighbor(pimpl->getBox(i + 1, j - 1, k - 1));
                                        me->addNeighbor(pimpl->getBox(i + 1, j - 1, k + 0));
                                        me->addNeighbor(pimpl->getBox(i + 1, j - 1, k + 1));
                                        me->addNeighbor(pimpl->getBox(i + 1, j + 0, k - 1));
                                        me->addNeighbor(pimpl->getBox(i + 1, j + 0, k + 0));
                                        me->addNeighbor(pimpl->getBox(i + 1, j + 0, k + 1));
                                        me->addNeighbor(pimpl->getBox(i + 1, j + 1, k - 1));
                                        me->addNeighbor(pimpl->getBox(i + 1, j + 1, k + 0));
                                        me->addNeighbor(pimpl->getBox(i + 1, j + 1, k + 1));
                                    }
                                }
                            }
                        }
                    }

                    if (pimpl->maxCutoff > 0) {

                        for (auto &&box : pimpl->boxes) {
                            box.particleIndices.clear();
                        }

                        auto it_pos = data.cbegin_positions();
                        long idx = 0;
                        const auto shift = readdy::model::Vec3(.5 * simBoxSize[0], .5 * simBoxSize[1],
                                                               .5 * simBoxSize[2]);
                        while (it_pos != data.cend_positions()) {
                            const auto pos_shifted = *it_pos + shift;
                            const long i = (const long) floor(pos_shifted[0] / pimpl->boxSize[0]);
                            const long j = (const long) floor(pos_shifted[1] / pimpl->boxSize[1]);
                            const long k = (const long) floor(pos_shifted[2] / pimpl->boxSize[2]);
                            auto box = pimpl->getBox(i, j, k);
                            if (box) {
                                box->particleIndices.push_back(idx);
                            }
                            ++idx;
                            ++it_pos;
                        }

                        for (auto &&box : pimpl->boxes) {
                            for (long i = 0; i < box.particleIndices.size(); ++i) {
                                const auto pI = box.particleIndices[i];
                                for (long j = i + 1; j < box.particleIndices.size(); ++j) {
                                    pimpl->pairs.emplace(pI, box.particleIndices[j]);
                                }

                                for (auto &&neighboringBox : box.neighboringBoxes) {
                                    for (const auto &pJ : neighboringBox->particleIndices) {
                                        pimpl->pairs.emplace(pI, pJ);
                                    }
                                }
                            }
                        }
                    }

                }

                iter_type NotThatNaiveSingleCPUNeighborList::begin() {
                    return pimpl->pairs.begin();
                }

                const_iter_type NotThatNaiveSingleCPUNeighborList::begin() const {
                    return cbegin();
                }

                const_iter_type NotThatNaiveSingleCPUNeighborList::cbegin() const {
                    return pimpl->pairs.cbegin();
                }

                iter_type NotThatNaiveSingleCPUNeighborList::end() {
                    return pimpl->pairs.end();
                }

                const_iter_type NotThatNaiveSingleCPUNeighborList::end() const {
                    return cend();
                }

                const_iter_type NotThatNaiveSingleCPUNeighborList::cend() const {
                    return pimpl->pairs.cend();
                }

                NotThatNaiveSingleCPUNeighborList::~NotThatNaiveSingleCPUNeighborList() = default;


            }
        }
    }
}
