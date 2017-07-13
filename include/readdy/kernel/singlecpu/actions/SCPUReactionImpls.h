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

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
namespace reactions {

class SCPUUncontrolledApproximation : public readdy::model::actions::reactions::UncontrolledApproximation {

public:
    SCPUUncontrolledApproximation(SCPUKernel *const kernel, scalar timeStep);

    virtual void perform() override;

    virtual void registerReactionScheme_11(const std::string &reactionName, reaction_11 fun) override;

    virtual void registerReactionScheme_12(const std::string &reactionName, reaction_12 fun) override;

    virtual void registerReactionScheme_21(const std::string &reactionName, reaction_21 fun) override;

    virtual void registerReactionScheme_22(const std::string &reactionName, reaction_22 fun) override;

protected:
    SCPUKernel *const kernel;
};

struct Event {
    using index_type = model::SCPUParticleData::index_t;
    using reaction_index_type = std::size_t;
    using particletype_t = readdy::model::Particle::type_type;
    unsigned int nEducts;
    unsigned int nProducts;
    index_type idx1, idx2;
    reaction_index_type reactionIdx;
    particletype_t t1, t2;
    scalar reactionRate;
    scalar cumulativeRate;

    Event(unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, scalar reactionRate,
          scalar cumulativeRate, reaction_index_type reactionIdx, particletype_t t1, particletype_t t2);

    friend std::ostream &operator<<(std::ostream &, const Event &);

};

class SCPUGillespie : public readdy::model::actions::reactions::Gillespie {
    using reaction_idx_t = Event::index_type;
public:

    SCPUGillespie(SCPUKernel *const kernel, scalar timeStep)
            : readdy::model::actions::reactions::Gillespie(timeStep), kernel(kernel) {};

    virtual void perform() override;

protected:
    SCPUKernel *const kernel;
};
}
}
}
}
}
