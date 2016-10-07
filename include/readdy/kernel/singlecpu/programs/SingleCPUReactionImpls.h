/**
 * << detailed description >>
 *
 * @file SingleCPUDefaultReactionProgram.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#ifndef READDY_MAIN_SINGLECPUDEFAULTREACTIONPROGRAM_H
#define READDY_MAIN_SINGLECPUDEFAULTREACTIONPROGRAM_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace programs {
namespace reactions {
class UncontrolledApproximation : public readdy::model::programs::reactions::UncontrolledApproximation {

public:
    UncontrolledApproximation(SingleCPUKernel const *const kernel);

    virtual void execute() override;

    virtual void registerReactionScheme_11(const std::string &name, reaction_11 fun) override;

    virtual void registerReactionScheme_12(const std::string &name, reaction_12 fun) override;

    virtual void registerReactionScheme_21(const std::string &name, reaction_21 fun) override;

    virtual void registerReactionScheme_22(const std::string &name, reaction_22 fun) override;

protected:
    SingleCPUKernel const *const kernel;
    std::unordered_map<short, reaction_11> mapping_11{};
    std::unordered_map<short, reaction_12> mapping_12{};
    std::unordered_map<short, reaction_21> mapping_21{};
    std::unordered_map<short, reaction_22> mapping_22{};
};

struct ReactionEvent {
    using index_type = std::size_t;
    unsigned int nEducts;
    unsigned int nProducts;
    index_type idx1, idx2;
    index_type reactionIdx;
    unsigned int t1, t2;
    double reactionRate;
    double cumulativeRate;

    ReactionEvent(unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, double reactionRate,
                  double cumulativeRate, index_type reactionIdx, unsigned int t1, unsigned int t2);

    friend std::ostream &operator<<(std::ostream &, const ReactionEvent &);

};

class Gillespie : public readdy::model::programs::reactions::Gillespie {
    using reaction_idx_t = ReactionEvent::index_type;
public:

    Gillespie(SingleCPUKernel const *const kernel) : kernel(kernel) {};

    virtual void execute() override {
        const auto &ctx = kernel->getKernelContext();
        auto data = kernel->getKernelStateModel().getParticleData();
        const auto &dist = ctx.getDistSquaredFun();
        const auto &fixPos = ctx.getFixPositionFun();

        double alpha = 0.0;
        auto events = gatherEvents(alpha);
        auto newParticles = handleEvents(std::move(events), alpha);

        // reposition particles to respect the periodic b.c.
        std::for_each(newParticles.begin(), newParticles.end(),
                      [&fixPos](readdy::model::Particle &p) { fixPos(p.getPos()); });

        // update data structure
        data->deactivateMarked();
        data->addParticles(newParticles);
    }

protected:
    virtual std::vector<ReactionEvent> gatherEvents(double &alpha);

    virtual std::vector<readdy::model::Particle> handleEvents(std::vector<ReactionEvent> events, double alpha);

    SingleCPUKernel const *const kernel;
};
}
}
}
}
}
#endif //READDY_MAIN_SINGLECPUDEFAULTREACTIONPROGRAM_H
