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
#include <readdy/kernel/cpu/util/scoped_thread.h>
#include <readdy/kernel/cpu/util/Config.h>
#include "ParticleData.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {
struct Neighbor {
    using index_t = readdy::model::Particle::id_type;
    index_t idx;
    double d2;

    Neighbor(const index_t idx, const double d2);
};

class NeighborList {

public:
    using box_index = unsigned short;
    using signed_box_index = typename std::make_signed<box_index>::type;
    using particle_index = Neighbor::index_t;
    using neighbor_t = Neighbor;
    using container_t = std::unordered_map<particle_index, std::vector<neighbor_t>>;
    using data_t = readdy::kernel::cpu::model::ParticleData;

    container_t pairs {};

    NeighborList(const readdy::model::KernelContext *const context, util::Config const *const config);

    virtual ~NeighborList();

    virtual void setupBoxes();

    virtual void setupNeighboringBoxes(const signed_box_index i, const signed_box_index j, const signed_box_index k);

    void clear();

    virtual void fillBoxes(const data_t &data);

    virtual void create(const data_t &data);

    void remove(const particle_index);

    void insert(const data_t &data, const particle_index);

protected:
    struct Box;
    const readdy::model::KernelContext *const ctx;

    using box_size_t = decltype(ctx->getBoxSize());

    box_size_t simBoxSize;

    std::vector<Box> boxes;
    std::array<box_index, 3> nBoxes{{0, 0, 0}};
    readdy::model::Vec3 boxSize{0, 0, 0};
    double maxCutoff = 0;
    util::Config const *const config;

    Box *getBox(const readdy::model::Particle::pos_type &pos);

    Box *getBox(signed_box_index i, signed_box_index j, signed_box_index k);

};
}
}
}
}
#endif //READDY_CPUKERNEL_NEIGHBORLIST_H
