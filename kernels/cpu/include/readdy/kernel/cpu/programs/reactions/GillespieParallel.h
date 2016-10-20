/**
 * << detailed description >>
 *
 * @file GillespieParallel.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_CPUKERNEL_GILLESPIEPARALLEL_H
#define READDY_CPUKERNEL_GILLESPIEPARALLEL_H


#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

class GillespieParallel : public readdy::model::programs::reactions::GillespieParallel {
    using kernel_t = readdy::kernel::cpu::CPUKernel;
    using vec_t = readdy::model::Vec3;
    using data_t = decltype(std::declval<kernel_t>().getKernelStateModel().getParticleData());
    using nl_t = const decltype(std::declval<kernel_t>().getKernelStateModel().getNeighborList());
    using ctx_t = std::remove_const<decltype(std::declval<kernel_t>().getKernelContext())>::type;
    using event_t = readdy::kernel::singlecpu::programs::reactions::ReactionEvent;
    using index_t = event_t::index_type;
    using particle_t = readdy::model::Particle;
public:
    GillespieParallel(kernel_t const *const kernel);

    ~GillespieParallel();

    virtual void execute() override;

    void clear();

    double getMaxReactionRadius() const;

    double getBoxWidth() const;

    unsigned int getLongestAxis() const;

    unsigned int getOtherAxis1() const;

    unsigned int getOtherAxis2() const;

    void setFilterEventsInAdvance(bool filterEventsInAdvance);

private:
    kernel_t const *const kernel;
    double maxReactionRadius = 0.0;
    double boxWidth = 0.0;
    unsigned int longestAxis;
    unsigned int otherAxis1, otherAxis2;
    bool filterEventsInAdvance = false;
    struct SlicedBox;
    std::vector<SlicedBox> boxes;
    bool approximateRate = true;

    /**
     * look for the longest axis and divide it into n_threads parts, yielding halo boxes.
     */
    void setupBoxes();

public:
    void setApproximateRate(bool approximateRate);

private:

    /**
     * Sorts particles into boxes and adds them to the problematic particles should they fall into
     * a halo region
     */
    void fillBoxes();

    /**
     * Executes the gillespie algorithm for each box and gives an update on problematic particles
     */
    void handleBoxReactions();

    template<typename ParticleCollection>
    void gatherEvents(const ParticleCollection &particles, const nl_t nl, const data_t data, double &alpha,
                      std::vector<GillespieParallel::event_t> &events) const;

    void findProblematicParticles(const unsigned long idx, const SlicedBox &box, ctx_t ctx,
                                  data_t data, nl_t nl, std::set<unsigned long> &problematic) const;

};
}

}
}
}
}

#endif //READDY_CPUKERNEL_GILLESPIEPARALLEL_H
