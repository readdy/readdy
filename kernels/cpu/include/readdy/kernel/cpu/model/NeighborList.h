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
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>
#include <readdy/kernel/cpu/util/ConfigUtils.h>
#include <readdy/kernel/cpu/util/ScopedThread.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace model {
                class NeighborList {
                public:
                    using container_t = std::unordered_map<unsigned long, std::vector<unsigned long>>;
                    NeighborList(const readdy::model::KernelContext *const context) : ctx(context) { }

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
                                for (unsigned long i = 0; i < nBoxes[0]; ++i) {
                                    for (unsigned long j = 0; j < nBoxes[1]; ++j) {
                                        for (unsigned long k = 0; k < nBoxes[2]; ++k) {
                                            setupNeighboringBoxes(i,j,k);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    virtual void setupNeighboringBoxes(unsigned long i, unsigned long j, unsigned long k) {
                        auto me = getBox(i, j, k);
                        for(int _i = 0; _i < 3; ++_i) {
                            for(int _j = 0; _j < 3; ++_j) {
                                for(int _k = 0; _k < 3; ++_k) {
                                    me->addNeighbor(getBox(_i, _j, _k));
                                }
                            }
                        }
                    }

                    virtual void fillBoxes(const singlecpu::model::SingleCPUParticleData &data) {
                        const auto simBoxSize = ctx->getBoxSize();
                        if (maxCutoff > 0) {

                            for (auto &&box : boxes) {
                                box.particleIndices.clear();
                            }
                            pairs->clear();

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
                                (*pairs)[idx] = std::vector<unsigned long>();
                                ++idx;
                                ++it_pos;
                            }


                            {
                                const auto size = boxes.size();
                                const std::size_t grainSize = size / util::getNThreads();

                                auto worker = [this](unsigned long begin, unsigned long end) {
                                    for(auto _b = begin; _b < end; ++_b) {
                                        const auto& box = boxes[_b];
                                        for (long i = 0; i < box.particleIndices.size(); ++i) {
                                            const auto pI = box.particleIndices[i];
                                            for (long j = 0; j < box.particleIndices.size(); ++j) {
                                                if(i != j) (*pairs)[pI].push_back(box.particleIndices[j]);
                                            }

                                            for (auto &&neighboringBox : box.neighboringBoxes) {
                                                for (const auto &pJ : neighboringBox->particleIndices) {
                                                    (*pairs)[pI].push_back(pJ);
                                                }
                                            }
                                        }
                                    }
                                };

                                std::vector<util::ScopedThread> threads;
                                threads.reserve(util::getNThreads());

                                for(auto i = 0; i < util::getNThreads()-1; ++i) {
                                    threads.push_back(util::ScopedThread(std::thread(worker, i*grainSize, (i+1)*grainSize)));
                                }
                                threads.push_back(util::ScopedThread(std::thread(worker, (util::getNThreads()-1)*grainSize, boxes.size())));
                            }

                                /*for (auto &&box : boxes) {
                                    for (long i = 0; i < box.particleIndices.size(); ++i) {
                                        const auto pI = box.particleIndices[i];
                                        for (long j = 0; j < box.particleIndices.size(); ++j) {
                                            if(i != j) (*pairs)[pI].push_back(box.particleIndices[j]);
                                        }

                                        for (auto &&neighboringBox : box.neighboringBoxes) {
                                            for (const auto &pJ : neighboringBox->particleIndices) {
                                                (*pairs)[pI].push_back(pJ);
                                            }
                                        }
                                    }
                                }*/


                        }
                    }

                    virtual void create(const readdy::kernel::singlecpu::model::SingleCPUParticleData &data) {
                        setupBoxes();
                        fillBoxes(data);
                    }

                    std::unique_ptr<container_t> pairs = std::make_unique<container_t>();

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
                    const readdy::model::KernelContext *const ctx;
                };
            }
        }
    }
}
#endif //READDY_MAIN_NEIGHBORLIST_H
