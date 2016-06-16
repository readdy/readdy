/**
 * << detailed description >>
 *
 * @file SingleCPUNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

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

                NaiveSingleCPUNeighborList &NaiveSingleCPUNeighborList::operator=(NaiveSingleCPUNeighborList &&rhs) = default;


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
                };

                NaiveSingleCPUNeighborList::~NaiveSingleCPUNeighborList() = default;

                NaiveSingleCPUNeighborList::NaiveSingleCPUNeighborList(NaiveSingleCPUNeighborList &&rhs) = default;


            }
        }
    }
}
