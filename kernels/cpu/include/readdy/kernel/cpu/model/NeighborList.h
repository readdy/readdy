/**
 * << detailed description >>
 *
 * @file NeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.07.16
 */

#ifndef READDY_CPUKERNEL_NEIGHBORLIST_H
#define READDY_CPUKERNEL_NEIGHBORLIST_H

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

    using index_t = readdy::model::Particle::id_type;

    const index_t idx;
    const double d2;

    NeighborListElement(const index_t idx, const double d2);
};

class NeighborList {
    struct Box;

public:
    using box_index = unsigned short;
    using signed_box_index = typename std::make_signed<box_index>::type;
    using neighbor_t = NeighborListElement;
    using container_t = std::unordered_map<NeighborListElement::index_t, std::vector<neighbor_t>>;

    std::unique_ptr<container_t> pairs = std::make_unique<container_t>();

    NeighborList(const readdy::model::KernelContext *const context, util::Config const *const config);

    virtual ~NeighborList();

    virtual void setupBoxes();

    virtual void setupNeighboringBoxes(const signed_box_index i, const signed_box_index j, const signed_box_index k);

    void clear();

    virtual void fillBoxes(const singlecpu::model::SingleCPUParticleData &data);

    virtual void create(const readdy::kernel::singlecpu::model::SingleCPUParticleData &data);

protected:

    std::vector<Box> boxes;
    std::array<box_index, 3> nBoxes{{0, 0, 0}};
    readdy::model::Vec3 boxSize{0, 0, 0};
    double maxCutoff = 0;
    util::Config const *const config;

    Box *getBox(signed_box_index i, signed_box_index j, signed_box_index k);

    const readdy::model::KernelContext *const ctx;
};
}
}
}
}
#endif //READDY_CPUKERNEL_NEIGHBORLIST_H
