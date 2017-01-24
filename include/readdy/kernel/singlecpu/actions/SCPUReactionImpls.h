/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file SingleCPUDefaultReactionProgram.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#ifndef READDY_MAIN_SINGLECPUDEFAULTREACTIONACTION_H
#define READDY_MAIN_SINGLECPUDEFAULTREACTIONACTION_H

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
namespace reactions {
class SCPUUncontrolledApproximation : public readdy::model::actions::reactions::UncontrolledApproximation {

public:
    SCPUUncontrolledApproximation(SCPUKernel const *const kernel, double timeStep);

    virtual void perform() override;

    virtual void registerReactionScheme_11(const std::string &name, reaction_11 fun) override;

    virtual void registerReactionScheme_12(const std::string &name, reaction_12 fun) override;

    virtual void registerReactionScheme_21(const std::string &name, reaction_21 fun) override;

    virtual void registerReactionScheme_22(const std::string &name, reaction_22 fun) override;

protected:
    SCPUKernel const *const kernel;
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

class SCPUGillespie : public readdy::model::actions::reactions::Gillespie {
    using reaction_idx_t = ReactionEvent::index_type;
public:

    SCPUGillespie(SCPUKernel const *const kernel, double timeStep)
            : readdy::model::actions::reactions::Gillespie(timeStep), kernel(kernel) {};

    virtual void perform() override {
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

    SCPUKernel const *const kernel;
};
}
}
}
}
}
#endif //READDY_MAIN_SINGLECPUDEFAULTREACTIONACTION_H
