/**
 * << detailed description >>
 *
 * @file SingleCPUNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace model {

                void NaiveSingleCPUNeighborList::create(const SingleCPUParticleData &data) {
                    pairs->clear();
                    for (size_t i = 0; i < data.size(); ++i) {
                        for (size_t j = i + 1; j < data.size(); ++j) {
                            pairs->emplace(i, j);
                        }
                    }
                }

                NaiveSingleCPUNeighborList &NaiveSingleCPUNeighborList::operator=(
                        NaiveSingleCPUNeighborList &&rhs) = default;


                NaiveSingleCPUNeighborList::NaiveSingleCPUNeighborList() = default;

                NaiveSingleCPUNeighborList::~NaiveSingleCPUNeighborList() = default;

                NaiveSingleCPUNeighborList::NaiveSingleCPUNeighborList(NaiveSingleCPUNeighborList &&rhs) = default;

                template<typename T>
                NotThatNaiveSingleCPUNeighborList<T>::Box::Box(long i, long j, long k, long id) : i(i), j(j), k(k), id(id) {
                };

                template<typename T>
                void NotThatNaiveSingleCPUNeighborList<T>::Box::addNeighbor(Box *box) {
                    if (box && box->id != id) neighboringBoxes.push_back(box);
                }

                template<typename T>
                typename NotThatNaiveSingleCPUNeighborList<T>::Box *NotThatNaiveSingleCPUNeighborList<T>::getBox(
                        long i, long j, long k
                ) {
                        const auto &periodic = ctx->getPeriodicBoundary();
                        if (periodic[0]) i = positive_modulo(i, nBoxes[0]);
                        else if (i < 0 || i >= nBoxes[0]) return nullptr;
                        if (periodic[1]) j = positive_modulo(j, nBoxes[1]);
                        else if (j < 0 || j >= nBoxes[1]) return nullptr;
                        if (periodic[2]) k = positive_modulo(k, nBoxes[2]);
                        else if (k < 0 || k >= nBoxes[2]) return nullptr;
                        return &boxes[k + j * nBoxes[2] + i * nBoxes[2] * nBoxes[1]];
                    }

                template<typename T>
                NotThatNaiveSingleCPUNeighborList<T>::NotThatNaiveSingleCPUNeighborList(
                        const readdy::model::KernelContext *const ctx) : ctx(ctx) {
                }

                template<typename T>
                void NotThatNaiveSingleCPUNeighborList<T>::create(const SingleCPUParticleData &data) {
                    setupBoxes();
                    fillBoxes(data);
                }

                template<typename T>
                void NotThatNaiveSingleCPUNeighborList<T>::setupBoxes() {
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
                        NotThatNaiveSingleCPUNeighborList<T>::maxCutoff = maxCutoff;
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

                template<typename T>
                void NotThatNaiveSingleCPUNeighborList<T>::fillBoxes(const SingleCPUParticleData &data) {
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
                                    super::pairs->emplace(pI, box.particleIndices[j]);
                                }

                                for (auto &&neighboringBox : box.neighboringBoxes) {
                                    for (const auto &pJ : neighboringBox->particleIndices) {
                                        super::pairs->emplace(pI, pJ);
                                    }
                                }
                            }
                        }
                    }
                }

                template<typename T>
                long NotThatNaiveSingleCPUNeighborList<T>::positive_modulo(long i, long n) const {
                    return (i % n + n) % n;
                }


                template<typename T>
                NotThatNaiveSingleCPUNeighborList<T>::~NotThatNaiveSingleCPUNeighborList() = default;

                template struct NotThatNaiveSingleCPUNeighborList<>;
            }
        }
    }
}
