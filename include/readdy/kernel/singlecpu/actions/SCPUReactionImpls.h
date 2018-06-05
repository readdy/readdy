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
#include "SCPUReactionUtils.h"

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
namespace reactions {

class SCPUUncontrolledApproximation : public readdy::model::actions::reactions::UncontrolledApproximation {

public:
    SCPUUncontrolledApproximation(SCPUKernel* kernel, scalar timeStep);

    void perform(const util::PerformanceNode &node) override;

    void registerReactionScheme_11(const std::string &reactionName, reaction_11 fun) override;

    void registerReactionScheme_12(const std::string &reactionName, reaction_12 fun) override;

    void registerReactionScheme_21(const std::string &reactionName, reaction_21 fun) override;

    void registerReactionScheme_22(const std::string &reactionName, reaction_22 fun) override;

protected:
    SCPUKernel *const kernel;
};

struct Event {
    using index_type = model::SCPUParticleData::entry_index;
    using reaction_index_type = std::size_t;
    std::uint8_t nEducts;
    std::uint8_t nProducts;
    index_type idx1, idx2;
    reaction_index_type reactionIndex;
    particle_type_type t1, t2;
    scalar rate;
    scalar cumulativeRate;

    Event(unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, scalar reactionRate,
          scalar cumulativeRate, reaction_index_type reactionIdx, particle_type_type t1, particle_type_type t2);

    friend std::ostream &operator<<(std::ostream &, const Event &);

};

class SCPUGillespie : public readdy::model::actions::reactions::Gillespie {
    using reaction_index = Event::index_type;
public:

    SCPUGillespie(SCPUKernel *const kernel, scalar timeStep)
            : readdy::model::actions::reactions::Gillespie(timeStep), kernel(kernel) {};

    void perform(const util::PerformanceNode &node) override;

protected:
    SCPUKernel *const kernel;
};

struct ParticleBackup {
    std::uint8_t nParticles; // either 1 or 2
    using index_type = model::SCPUParticleData::entry_index;
    index_type idx1, idx2;
    particle_type_type t1, t2;
    Vec3 pos1, pos2;

    ParticleBackup(Event event, const readdy::model::actions::reactions::ReversibleReactionConfig *revReaction,
                   const readdy::model::reactions::Reaction *reaction, const scpu_data *data);
};

class SCPUDetailedBalance : public readdy::model::actions::reactions::DetailedBalance {
    using scpu_data = readdy::kernel::scpu::model::SCPUParticleData;
    using fix_pos = readdy::model::Context::fix_pos_fun;
    using reaction_type = readdy::model::reactions::ReactionType;
public:
    SCPUDetailedBalance(SCPUKernel *const kernel, scalar timeStep)
            : readdy::model::actions::reactions::DetailedBalance(timeStep), kernel(kernel) {
        kernel->context().configure();
        searchReversibleReactions(kernel->context());
    };

    void perform(const util::PerformanceNode &node) override;

protected:
    SCPUKernel *const kernel;

    // calculate first-order interactions and second-order non-bonded interactions
    void calculateEnergies(const util::PerformanceNode &node);

    std::pair<model::SCPUParticleData::entries_update, scalar>
    performReversibleReactionEvent(const Event &event,
                                   const readdy::model::actions::reactions::ReversibleReactionConfig *reversibleReaction,
                                   const readdy::model::reactions::Reaction *reaction, reaction_record *record);

    model::SCPUParticleData::entries_update
    generateBackwardUpdate(const ParticleBackup &particleBackup,
                           const std::vector<model::SCPUParticleData::entry_index> &updateRecord) const;

    std::pair<const readdy::model::actions::reactions::ReversibleReactionConfig *, const readdy::model::reactions::Reaction *>
    findReversibleReaction(const Event &event);
};

}
}
}
}
}
