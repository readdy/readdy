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
                    ParticleIndexPair(size_t idx1, size_t idx2){
                        if(idx1 < idx2) {
                            ParticleIndexPair::idx1 = idx1;
                            ParticleIndexPair::idx2 = idx2;
                        } else if(idx1 > idx2){
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

                    friend std::ostream& operator<<(std::ostream& os, const ParticleIndexPair &pip) {
                        os << "ParticleIndexPair(" << pip.idx1 << ", " << pip.idx2 <<")";
                        return os;
                    }
                };

                struct ParticleIndexPairHasher {
                    size_t operator()(const ParticleIndexPair &pip) const {
                        return hash_value(pip);
                    }
                };

                using iter_type = std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher>::iterator;
                using const_iter_type = std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher>::const_iterator;

                struct SingleCPUNeighborList {
                    virtual void create(const SingleCPUParticleData &data) = 0;

                    virtual iter_type begin() = 0;
                    virtual const_iter_type begin() const = 0;
                    virtual const_iter_type cbegin() const = 0;
                    virtual iter_type end() = 0;
                    virtual const_iter_type end() const = 0;
                    virtual const_iter_type cend() const = 0;
                };

                struct NaiveSingleCPUNeighborList : public SingleCPUNeighborList{
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

                    virtual iter_type begin() override;

                    virtual const_iter_type begin() const override;

                    virtual const_iter_type cbegin() const override;

                    virtual iter_type end() override;

                    virtual const_iter_type end() const override;

                    virtual const_iter_type cend() const override;


                protected:
                    struct Impl;
                    std::unique_ptr<Impl> pimpl;
                };

                struct NotThatNaiveSingleCPUNeighborList : public SingleCPUNeighborList {

                    NotThatNaiveSingleCPUNeighborList(const readdy::model::KernelContext *const ctx);
                    virtual ~NotThatNaiveSingleCPUNeighborList();

                    virtual void create(const SingleCPUParticleData &data) override;
                    virtual void setupBoxes();
                    virtual void fillBoxes(const SingleCPUParticleData &data);

                    virtual iter_type begin() override;
                    virtual const_iter_type begin() const override;
                    virtual const_iter_type cbegin() const override;
                    virtual iter_type end() override;
                    virtual const_iter_type end() const override;
                    virtual const_iter_type cend() const override;

                protected:
                    struct Box {
                        std::vector<Box *> neighboringBoxes{};
                        std::vector<long> particleIndices{};
                        long i, j, k;
                        long id = 0;

                        Box(long i, long j, long k, long id);
                        void addNeighbor(Box *box);
                    };
                    const readdy::model::KernelContext *ctx;
                    std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher> pairs{};
                    std::vector<Box> boxes{};
                    std::array<int, 3> nBoxes{{0, 0, 0}};
                    readdy::model::Vec3 boxSize{0, 0, 0};
                    double maxCutoff = 0;

                    long positive_modulo(long i, long n) const;
                    Box *getBox(long i, long j, long k);

                };

            }
        }
    }
}

#endif //READDY_MAIN_SINGLECPUNEIGHBORLIST_H
