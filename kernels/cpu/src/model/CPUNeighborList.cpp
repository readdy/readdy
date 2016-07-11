/**
 * << detailed description >>
 *
 * @file CPUNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */

#include <readdy/kernel/cpu/model/CPUNeighborList.h>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace model {
                NotThatNaiveCPUNeighborList::NotThatNaiveCPUNeighborList(const readdy::model::KernelContext *const ctx)
                        : readdy::kernel::singlecpu::model::NotThatNaiveSingleCPUNeighborList(ctx) { }

                void NotThatNaiveCPUNeighborList::fillBoxes(const readdy::kernel::singlecpu::model::SingleCPUParticleData &data) {
                    const auto simBoxSize = ctx->getBoxSize();
                    if (maxCutoff > 0) {

                        for (auto &&box : boxes) {
                            box.particleIndices.clear();
                        }

                        auto it_pos = data.cbegin_positions();
                        long idx = 0;
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
                                    pairs.emplace(pI, box.particleIndices[j]);
                                }

                                for (auto &&neighboringBox : box.neighboringBoxes) {
                                    for (const auto &pJ : neighboringBox->particleIndices) {
                                        pairs.emplace(pI, pJ);
                                    }
                                }
                            }
                        }
                    }
                }


            }
        }
    }
}