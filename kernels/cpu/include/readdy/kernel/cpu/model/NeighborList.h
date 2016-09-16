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
struct NeighborListElement {
    const unsigned long idx;
    const double d2;
    NeighborListElement(const unsigned long idx, const double d2);
};

class NeighborList {
    struct Box;
    using box_t = Box;

public:
    using neighbor_t = NeighborListElement;
    using container_t = std::unordered_map<unsigned long, std::vector<neighbor_t>>;

    std::unique_ptr<container_t> pairs = std::make_unique<container_t>();

    NeighborList(const readdy::model::KernelContext *const context, util::Config const *const config);
    virtual ~NeighborList();

    virtual void setupBoxes();

    virtual void setupNeighboringBoxes(unsigned long i, unsigned long j, unsigned long k);

    void clear();

    virtual void fillBoxes(const singlecpu::model::SingleCPUParticleData &data);

    virtual void create(const readdy::kernel::singlecpu::model::SingleCPUParticleData &data);

protected:

    std::vector<Box> boxes;
    std::array<unsigned int, 3> nBoxes{{0, 0, 0}};
    readdy::model::Vec3 boxSize{0, 0, 0};
    double maxCutoff = 0;
    util::Config const *const config;

    long positive_modulo(long i, long n) const;

    Box *getBox(long i, long j, long k);

    const readdy::model::KernelContext *const ctx;
};
}
}
}
}
#endif //READDY_MAIN_NEIGHBORLIST_H
