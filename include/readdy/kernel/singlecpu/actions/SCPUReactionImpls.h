/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * Reaction handler declaration specific to Single CPU kernel. Defines additional structs: a reaction Event
 * to be used by the reaction handlers, and ParticleBackup to be used by the DetailedBalance reaction handler.
 *
 * @file SCPUReactionImpls.h
 * @brief Single CPU kernel declaration of reaction handlers
 * @author clonker
 * @author chrisfroe
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

    void perform() override;

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
    ParticleTypeId t1, t2;
    scalar rate;
    scalar cumulativeRate;

    Event(unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, scalar reactionRate,
          scalar cumulativeRate, reaction_index_type reactionIdx, ParticleTypeId t1, ParticleTypeId t2);

    friend std::ostream &operator<<(std::ostream &, const Event &);

};

class SCPUGillespie : public readdy::model::actions::reactions::Gillespie {
    using reaction_index = Event::index_type;
public:

    SCPUGillespie(SCPUKernel *const kernel, scalar timeStep)
            : readdy::model::actions::reactions::Gillespie(timeStep), kernel(kernel) {};

    void perform() override;

protected:
    SCPUKernel *const kernel;
};

struct ParticleBackup {
    std::uint8_t nParticles; // either 1 or 2
    using index_type = model::SCPUParticleData::entry_index;
    index_type idx1, idx2;
    ParticleTypeId t1, t2;
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
        searchReversibleReactions(kernel->context());
    };

    void perform() override;

protected:
    SCPUKernel *const kernel;

    // calculate first-order interactions and second-order non-bonded interactions
    void calculateEnergies();

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
