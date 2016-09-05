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
#include <readdy/kernel/cpu/util/ScopedThread.h>
#include <readdy/kernel/cpu/util/Config.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace model {
                class NeighborList {
                    using box_t = readdy::kernel::singlecpu::model::Box;

                public:
                    using container_t = std::unordered_map<unsigned long, std::vector<unsigned long>>;
                    NeighborList(const readdy::model::KernelContext *const context, util::Config const*const config)
                            : ctx(context), config(config) { }

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
                        for(int _i = -1; _i < 2; ++_i) {
                            for(int _j = -1; _j < 2; ++_j) {
                                for(int _k = -1; _k < 2; ++_k) {
                                    // don't add me as neighbor to myself
                                    if(!(_i == 0 && _j == 0 && _k == 0)) {
                                        me->addNeighbor(getBox(i + _i, j + _j, k + _k));
                                    }
                                }
                            }
                        }
                    }

                    void clear() {
                        boxes.clear();
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
                                const std::size_t grainSize = size / config->nThreads;

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
                                threads.reserve(config->nThreads);

                                for(auto i = 0; i < config->nThreads-1; ++i) {
                                    threads.push_back(util::ScopedThread(std::thread(worker, i*grainSize, (i+1)*grainSize)));
                                }
                                threads.push_back(util::ScopedThread(std::thread(worker, (config->nThreads-1)*grainSize, boxes.size())));
                            }
                        }
                    }

                    virtual void create(const readdy::kernel::singlecpu::model::SingleCPUParticleData &data) {
                        setupBoxes();
                        fillBoxes(data);
                    }

                    std::unique_ptr<container_t> pairs = std::make_unique<container_t>();

                    const std::vector<box_t> &getBoxes() const {
                        return boxes;
                    }

                protected:

                    std::vector<box_t> boxes{};
                    std::array<int, 3> nBoxes{{0, 0, 0}};
                    readdy::model::Vec3 boxSize{0, 0, 0};
                    double maxCutoff = 0;
                    util::Config const*const config;

                    long positive_modulo(long i, long n) const {
                        return (i % n + n) % n;
                    }
                    box_t *getBox(long i, long j, long k) {
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
