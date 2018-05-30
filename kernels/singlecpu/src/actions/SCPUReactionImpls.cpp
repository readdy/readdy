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
 * @file SingleCPUDefaultReactionProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <readdy/kernel/singlecpu/actions/SCPUReactionImpls.h>
#include <readdy/kernel/singlecpu/actions/SCPUReactionUtils.h>
#include <readdy/common/algorithm.h>

#include <utility>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {

namespace reactions {

Event::Event(unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, scalar reactionRate,
             scalar cumulativeRate, reaction_index_type reactionIdx, particle_type_type t1, particle_type_type t2)
        : nEducts(nEducts), nProducts(nProducts), idx1(idx1), idx2(idx2), rate(reactionRate),
          cumulativeRate(cumulativeRate), reactionIndex(reactionIdx), t1(t1), t2(t2) {
}

std::ostream &operator<<(std::ostream &os, const Event &evt) {
    os << "Event(" << evt.idx1 << "[type=" << evt.t1 << "]";
    if (evt.nEducts == 2) {
        os << " + " << evt.idx2 << "[type=" << evt.t2 << "]";
    }
    os << ", rate=" << evt.rate << ", cumulativeRate=" << evt.cumulativeRate
       << ", reactionIdx=" << evt.reactionIndex;
    return os;
}

using event_t = Event;

SCPUUncontrolledApproximation::SCPUUncontrolledApproximation(SCPUKernel *const kernel, scalar timeStep)
        : readdy::model::actions::reactions::UncontrolledApproximation(timeStep), kernel(kernel) {
}

template<bool approximated>
bool performReactionEvent(const scalar rate, const scalar timeStep) {
    if (approximated) {
        return readdy::model::rnd::uniform_real() < rate * timeStep;
    } else {
        return readdy::model::rnd::uniform_real() < 1 - std::exp(-rate * timeStep);
    }
}


inline bool shouldPerformEvent(const scalar rate, const scalar timestep, bool approximated) {
    return approximated ? performReactionEvent<true>(rate, timestep) : performReactionEvent<false>(rate, timestep);
}

std::vector<event_t> findEvents(const SCPUKernel *const kernel, scalar dt, bool approximateRate = false) {
    std::vector<event_t> eventsUpdate;
    auto &stateModel = kernel->getSCPUKernelStateModel();
    auto &data = *stateModel.getParticleData();
    const auto &box = kernel->context().boxSize().data();
    const auto &pbc = kernel->context().periodicBoundaryConditions().data();
    {
        std::size_t idx = 0;
        for (const auto &e : data) {
            if (!e.is_deactivated()) {
                // order 1
                const auto &reactions = kernel->context().reactions().order1ByType(e.type);
                for (auto it_reactions = reactions.begin(); it_reactions != reactions.end(); ++it_reactions) {
                    const auto rate = (*it_reactions)->rate();
                    if (rate > 0 && shouldPerformEvent(rate, dt, approximateRate)) {
                        Event evt{1, (*it_reactions)->nProducts(), idx, idx, rate, 0,
                                  static_cast<std::size_t>(it_reactions - reactions.begin()),
                                  e.type, 0};
                        eventsUpdate.push_back(evt);
                    }
                }
            }
            ++idx;
        }
    }
    const auto &nl = *stateModel.getNeighborList();
    const auto &context = kernel->context();
    for(auto cell = 0_z; cell < nl.nCells(); ++cell) {
        for(auto it = nl.particlesBegin(cell); it != nl.particlesEnd(cell); ++it) {
            auto pidx = *it;
            const auto &entry = data.entry_at(*it);
            if(!entry.deactivated) {
                nl.forEachNeighbor(it, cell, [&](const std::size_t neighborIdx) {
                    const auto &neighbor = data.entry_at(neighborIdx);
                    if(!neighbor.deactivated) {
                        const auto &reactions = context.reactions().order2ByType(entry.type, neighbor.type);
                        if (!reactions.empty()) {
                            const auto distSquared = bcs::distSquared(neighbor.position(), entry.position(), box, pbc);
                            for (auto it_reactions = reactions.begin(); it_reactions < reactions.end(); ++it_reactions) {
                                const auto &react = *it_reactions;
                                const auto rate = react->rate();
                                if (rate > 0 && distSquared < react->eductDistanceSquared()
                                    && shouldPerformEvent(rate, dt, approximateRate)) {
                                    const auto reaction_index = static_cast<event_t::reaction_index_type>(it_reactions -
                                                                                                          reactions.begin());
                                    eventsUpdate.emplace_back(2, react->nProducts(), pidx, neighborIdx, rate, 0,
                                                              reaction_index, entry.type, neighbor.type);
                                }
                            }
                        }
                    }
                });
            }
        }
    }
    return eventsUpdate;
}

void SCPUUncontrolledApproximation::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    const auto &ctx = kernel->context();
    SCPUStateModel &stateModel = kernel->getSCPUKernelStateModel();
    if(ctx.recordReactionsWithPositions()) {
        stateModel.reactionRecords().clear();
    }
    if(ctx.recordReactionCounts()) {
        stateModel.resetReactionCounts();
    }
    auto &data = *stateModel.getParticleData();
    auto events = findEvents(kernel, timeStep, false);

    // shuffle reactions
    std::shuffle(events.begin(), events.end(), std::mt19937(std::random_device()()));

    // execute reactions
    {
        scpu_data::new_entries newParticles{};
        std::vector<scpu_data::entry_index> decayedEntries{};

        readdy::model::reactions::Reaction *reaction;

        for (auto it = events.begin(); it != events.end(); ++it) {
            auto &event = *it;
            if (event.cumulativeRate == 0) {
                auto entry1 = event.idx1;
                if (event.nEducts == 1) {
                    reaction = ctx.reactions().order1ByType(event.t1)[event.reactionIndex];
                    if(ctx.recordReactionsWithPositions()) {
                        reaction_record record;
                        record.id = reaction->id();
                        performReaction(data, entry1, entry1, newParticles, decayedEntries, reaction, ctx, &record);
                        stateModel.reactionRecords().push_back(record);
                    } else {
                        performReaction(data, entry1, entry1, newParticles, decayedEntries, reaction, ctx, nullptr);
                    }
                    for (auto _it2 = it + 1; _it2 != events.end(); ++_it2) {
                        if (_it2->idx1 == entry1 || _it2->idx2 == entry1) {
                            _it2->cumulativeRate = 1;
                        }
                    }
                } else {
                    reaction = ctx.reactions().order2ByType(event.t1, event.t2)[event.reactionIndex];
                    if(ctx.recordReactionsWithPositions()) {
                        reaction_record record;
                        record.id = reaction->id();
                        performReaction(data, entry1, event.idx2, newParticles, decayedEntries, reaction, ctx, &record);
                        stateModel.reactionRecords().push_back(record);
                    } else {
                        performReaction(data, entry1, event.idx2, newParticles, decayedEntries, reaction, ctx, nullptr);
                    }
                    for (auto _it2 = it + 1; _it2 != events.end(); ++_it2) {
                        if (_it2->idx1 == entry1 || _it2->idx2 == entry1 ||
                            _it2->idx1 == event.idx2 || _it2->idx2 == event.idx2) {
                            _it2->cumulativeRate = 1;
                        }
                    }
                }
                // update reaction counts if necessary
                if(ctx.recordReactionCounts()) {
                    stateModel.reactionCounts().at(reaction->id())++;
                }
            }
        }
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));
    }
}

void SCPUUncontrolledApproximation::registerReactionScheme_11(const std::string &reactionName,
                                                              readdy::model::actions::reactions::UncontrolledApproximation::reaction_11) {

}

void SCPUUncontrolledApproximation::registerReactionScheme_12(const std::string &reactionName,
                                                              readdy::model::actions::reactions::UncontrolledApproximation::reaction_12) {

}

void SCPUUncontrolledApproximation::registerReactionScheme_21(const std::string &reactionName,
                                                              readdy::model::actions::reactions::UncontrolledApproximation::reaction_21) {

}

void SCPUUncontrolledApproximation::registerReactionScheme_22(const std::string &reactionName,
                                                              readdy::model::actions::reactions::UncontrolledApproximation::reaction_22) {

}

void gatherEvents(SCPUKernel const *const kernel,
                  const readdy::kernel::scpu::model::CellLinkedList &nl, const scpu_data &data, scalar &alpha,
                  std::vector<event_t> &events) {
    {
        std::size_t index = 0;
        for (const auto &entry : data) {
            if (!entry.is_deactivated()) {
                const auto &reactions = kernel->context().reactions().order1ByType(entry.type);
                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    const auto rate = (*it)->rate();
                    if (rate > 0) {
                        alpha += rate;
                        events.emplace_back(1, (*it)->nProducts(), index, 0, rate, alpha,
                                            static_cast<event_t::reaction_index_type>(it - reactions.begin()),
                                            entry.type, 0);
                    }
                }
            }
            ++index;
        }
    }
    const auto box = kernel->context().boxSize().data();
    const auto pbc = kernel->context().periodicBoundaryConditions().data();
    for (auto cell = 0_z; cell < nl.nCells(); ++cell) {
        for (auto it = nl.particlesBegin(cell); it != nl.particlesEnd(cell); ++it) {
            auto &entry = data.entry_at(*it);
            if (!entry.deactivated) {
                nl.forEachNeighbor(it, cell, [&](const std::size_t nIdx) {
                    const auto &neighbor = data.entry_at(nIdx);
                    if (!neighbor.deactivated) {
                        const auto &reactions = kernel->context().reactions().order2ByType(entry.type,
                                                                                           neighbor.type);
                        if (!reactions.empty()) {
                            const auto distSquared = bcs::distSquared(neighbor.position(), entry.position(), box, pbc);
                            for (auto itR = reactions.begin(); itR < reactions.end(); ++itR) {
                                const auto &react = *itR;
                                const auto rate = react->rate();
                                if (rate > 0 && distSquared < react->eductDistanceSquared()) {
                                    alpha += rate;
                                    events.emplace_back(2, react->nProducts(), *it, nIdx,
                                                        rate, alpha,
                                                        static_cast<event_t::reaction_index_type>(
                                                                itR - reactions.begin()
                                                        ),
                                                        entry.type, neighbor.type);
                                }
                            }
                        }
                    }
                });
            }
        }
    }
}

scpu_data::entries_update handleEventsGillespie(
        SCPUKernel *const kernel, scalar timeStep,
        bool filterEventsInAdvance, bool approximateRate,
        std::vector<event_t> &&events) {
    scpu_data::new_entries newParticles{};
    std::vector<scpu_data::entry_index> decayedEntries{};

    if (!events.empty()) {
        const auto &ctx = kernel->context();
        auto &model = kernel->getSCPUKernelStateModel();
        const auto data = kernel->getSCPUKernelStateModel().getParticleData();
        /**
         * Handle gathered reaction events
         */
        {
            auto shouldEval = [&](const event_t &event) {
                return filterEventsInAdvance || shouldPerformEvent(event.rate, timeStep, approximateRate);
            };

            auto depending = [&](const event_t &e1, const event_t &e2) {
                return (e1.idx1 == e2.idx1 || (e2.nEducts == 2 && e1.idx1 == e2.idx2)
                        || (e1.nEducts == 2 && (e1.idx2 == e2.idx1 || (e2.nEducts == 2 && e1.idx2 == e2.idx2))));
            };

            auto eval = [&](const event_t &event) {
                auto entry1 = event.idx1;
                if (event.nEducts == 1) {
                    auto reaction = ctx.reactions().order1ByType(event.t1)[event.reactionIndex];
                    if(ctx.recordReactionsWithPositions()) {
                        reaction_record record;
                        record.id = reaction->id();
                        performReaction(*data, entry1, entry1, newParticles, decayedEntries, reaction, ctx,
                                        &record);
                        model.reactionRecords().push_back(record);
                    } else {
                        performReaction(*data, entry1, entry1, newParticles, decayedEntries, reaction, ctx,
                                        nullptr);
                    }
                    if(ctx.recordReactionCounts()) {
                        model.reactionCounts().at(reaction->id())++;
                    }
                } else {
                    auto reaction = ctx.reactions().order2ByType(event.t1, event.t2)[event.reactionIndex];
                    if(ctx.recordReactionsWithPositions()) {
                        reaction_record record;
                        record.id = reaction->id();
                        performReaction(*data, entry1, event.idx2, newParticles, decayedEntries, reaction, ctx, &record);
                        model.reactionRecords().push_back(record);
                    } else {
                        performReaction(*data, entry1, event.idx2, newParticles, decayedEntries, reaction, ctx, nullptr);
                    }
                    if(ctx.recordReactionCounts()) {
                        model.reactionCounts().at(reaction->id())++;
                    }
                }
            };

            algo::performEvents(events, shouldEval, depending, eval);

        }
    }
    return std::make_pair(std::move(newParticles), std::move(decayedEntries));
}

void SCPUGillespie::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    const auto &ctx = kernel->context();
    if(ctx.reactions().nOrder1() == 0 && ctx.reactions().nOrder2() == 0) {
        return;
    }
    auto &stateModel = kernel->getSCPUKernelStateModel();
    if(ctx.recordReactionsWithPositions()) stateModel.reactionRecords().clear();
    if(ctx.recordReactionCounts()) {
        stateModel.resetReactionCounts();
    }
    auto data = stateModel.getParticleData();
    const auto nl = stateModel.getNeighborList();

    scalar alpha = 0.0;
    std::vector<event_t> events;
    gatherEvents(kernel, *nl, *data, alpha, events);
    auto particlesUpdate = handleEventsGillespie(kernel, timeStep, false, false, std::move(events));

    // update data structure
    data->update(std::move(particlesUpdate));
}

ParticleBackup::ParticleBackup(unsigned int nParticles, ParticleBackup::index_type idx1,
                               ParticleBackup::index_type idx2, particle_type_type t1, particle_type_type t2, Vec3 pos1,
                               Vec3 pos2)
        : nParticles(nParticles), idx1(idx1), idx2(idx2), t1(t1), t2(t2), pos1(pos1), pos2(pos2) {
    if (!(nParticles == 1 || nParticles == 2)) {
        throw std::runtime_error(fmt::format("nParticles in ParticleBackup must be either 1 or 2, but was {}", nParticles));
    }
}

void SCPUDetailedBalance::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    log::trace("SCPUDetailedBalance::perform(){}{}", readdy::util::str::newline, readdy::util::str::newline);
    const auto &ctx = kernel->context();
    if (firstPerform) {
        searchReversibleReactions(ctx);
        firstPerform = false;
    }
    if(ctx.reactions().nOrder1() == 0 && ctx.reactions().nOrder2() == 0) {
        return;
    }
    auto &stateModel = kernel->getSCPUKernelStateModel();
    if(ctx.recordReactionCounts()) {
        stateModel.resetReactionCounts();
    }

    auto data = stateModel.getParticleData();
    const auto nl = stateModel.getNeighborList();

    scalar alpha = 0.0;
    std::vector<event_t> events;
    gatherEvents(kernel, *nl, *data, alpha, events);
    // @todo shuffle events

    calculateForcesEnergies();
    std::size_t nDeactivated = 0;
    const std::size_t nEvents = events.size();
    while (nDeactivated < nEvents) {
        log::trace("Handle next event, still {} events left", events.size());
        const auto event = events.back(); // copy on purpose, due to pop_back later in this loop
        if (not shouldPerformEvent(event.rate, timeStep, false)) {
            log::trace("Should not be performed, next one.");
            events.pop_back();
            nDeactivated++;
            continue;
        }
        bool isFusionReaction = (event.nEducts == 2);
        const auto &revReaction = _reversibleReactions.front(); // @todo search actual reversibleReactionConfig
        const auto particleBackup = ParticleBackup(event.nEducts, event.idx1, event.idx2, event.t1, event.t2, data->entry_at(event.idx1).position(), data->entry_at(event.idx2).position());
        auto energyBefore = stateModel.energy();
        model::SCPUParticleData::entries_update forwardUpdate;
        scalar interactionEnergy;
        std::tie(forwardUpdate, interactionEnergy) = performEvent(*data, event, ctx.recordReactionCounts());
        const auto updateRecord = data->update(std::move(forwardUpdate));
        auto backwardUpdate = generateBackwardUpdate(particleBackup, updateRecord);
        stateModel.updateNeighborList();
        calculateForcesEnergies();
        const auto deltaEnergy = (stateModel.energy() - interactionEnergy) - energyBefore;
        // get a-value from reversibleReactionConfig depending on direction of reaction
        const auto aValue = revReaction.acceptancePrefactor;
        auto acceptance = 1.;
        if (isFusionReaction) {
            acceptance = std::min(1., aValue * std::exp(-1. * deltaEnergy / ctx.kBT()));
        } else {
            acceptance = std::min(1., 1./aValue * std::exp(-1. * deltaEnergy / ctx.kBT()));
        }
        log::trace("Acceptance for current event is {}", acceptance);

        if (readdy::model::rnd::uniform_real() < acceptance) {
            // accept/do-nothing
            log::trace("accept!");
        } else {
            // reject/rollback
            log::trace("reject! subtract reaction count again");
            if (ctx.recordReactionCounts()) {
                if (isFusionReaction) {
                    const auto &reaction = ctx.reactions().order2ByType(event.t1, event.t2)[event.reactionIndex];
                    stateModel.reactionCounts().at(reaction->id())--;
                } else {
                    const auto &reaction = ctx.reactions().order1ByType(event.t1)[event.reactionIndex];
                    stateModel.reactionCounts().at(reaction->id())--;
                }
            }
            data->update(std::move(backwardUpdate));
            stateModel.updateNeighborList();
            calculateForcesEnergies();
            if (energyBefore != stateModel.energy()) {
                log::warn(
                        "reaction move was rejected but energy of state is different after rollback, was {}, is now {}",
                        energyBefore, stateModel.energy());
            }
        }
        events.pop_back();
        nDeactivated++;
        // remove all events involving idx1 or idx2
        const auto isRelevant = [&isFusionReaction](const event_t &testedEvent, const event_t &performedEvent) -> bool {
            if (isFusionReaction) {
                // there are two educts in performedEvent
                const auto decayedIdx1 = performedEvent.idx1;
                const auto decayedIdx2 = performedEvent.idx2;
                // check ANY occurrence of either decayedIdx1 or 2 in testedEvent
                if (testedEvent.nEducts == 1) {
                    // idx2 must not be checked, it might be set to anything
                    return (decayedIdx1 == testedEvent.idx1 || decayedIdx2 == testedEvent.idx1);
                } else if (testedEvent.nEducts == 2) {
                    return ((decayedIdx1 == testedEvent.idx1 || decayedIdx1 == testedEvent.idx2) or
                            (decayedIdx2 == testedEvent.idx1 || decayedIdx2 == testedEvent.idx2));
                } else {
                    throw std::runtime_error("unreachable");
                }

            } else {
                // there is only one educt in performedEvent
                const auto decayedIdx1 = performedEvent.idx1;
                if (testedEvent.nEducts == 1) {
                    // idx2 must not be checked, it might be set to anything
                    return (decayedIdx1 == testedEvent.idx1);
                } else if (testedEvent.nEducts == 2) {
                    return (decayedIdx1 == testedEvent.idx1 || decayedIdx1 == testedEvent.idx2);
                } else {
                    throw std::runtime_error("unreachable");
                }
            }
        };
        auto it = events.begin();
        auto end = events.end();
        size_t nRelevant = 0;
        while (it < end) {
            if (isRelevant(*it, event)) {
                ++nDeactivated;
                ++nRelevant;
                it = events.erase(it);
                end = events.end();
            } else {
                it++;
            }
        }
        log::trace("found and removed {} relevant events", nRelevant);
    }

    // Then calculate for each relevant fusion/fission pair:
    // acceptance forward, acceptance backward, array of radii and array of cumulatives to draw from

    // gather events, don't evaluate probabilities yet
    // for each event
        // energyBefore = stateModel.energy()
        // interactionEnergy is positive when positive energy was gained in fission
        // forwardUpdate, interactionEnergy = performEvent()
        // the record contains the indices of the new particles, len=2 for fission, len=1 for fusion
        // when rejecting, remove these and restore the oldstate from information from the event
        // updateRecord = data->update(forwardUpdate)
        // backwardUpdate = generateBackwardUpdate(forwardUpdate, updateRecord)
        // stateModel.updateNL()
        // stateModel.calculateForcesEnergies()
        // delta = (stateModel.energy() - interactionEnergy) - energyBefore
        // acceptance = min{1, a_value * e^{-beta delta}}
        // if uniform < acceptance
            // accept !
            // deactivate events whose educts have disappeared (including the just handled one)
            // therefore misuse the cumulativeRate as in Uncontrolled
        // else
            // reject ! apply backwardUpdate
            // data->update(backwardUpdate)
            // stateModel.updateNL()
            // stateModel.calculateForcesEnergies()
            // assert stateModel.energy() == energyBefore

}

model::SCPUParticleData::entries_update SCPUDetailedBalance::generateBackwardUpdate(const ParticleBackup &particleBackup, const std::vector<model::SCPUParticleData::entry_index> &updateRecord) const {
    std::vector<scpu_data::entry_index> decayedEntries = updateRecord; // the entries that were created by the forward update
    scpu_data::new_entries newParticles{};
    if (particleBackup.nParticles == 1) {
        // @fixme order of particles w.r.t. old indices might be mixed up,
        // but events referencing the same indices will have to be removed anyway, irrespective of the event
        // being accepted or rejected
        readdy::model::Particle p(particleBackup.pos1, particleBackup.t1);
        newParticles.emplace_back(p);
    } else if (particleBackup.nParticles == 2) {
        readdy::model::Particle p1(particleBackup.pos1, particleBackup.t1);
        readdy::model::Particle p2(particleBackup.pos2, particleBackup.t2);
        newParticles.emplace_back(p1);
        newParticles.emplace_back(p2);
    } else {
        throw std::runtime_error("Not gonna happen");
    }
    log::trace("generateBackwardUpdate: newParticles.size() {}, decayedEntries.size() {}", newParticles.size(),
               decayedEntries.size());
    log::trace("newParticles are ..");
    for (const auto &entry : newParticles) {
        log::trace("entry.position() {}, entry.type {}, entry.id {}, entry.deactivated {}", entry.position(),
                   entry.type, entry.id, entry.deactivated);
    }
    return std::make_pair(std::move(newParticles), decayedEntries);
}

std::pair<model::SCPUParticleData::entries_update, scalar> SCPUDetailedBalance::performEvent(scpu_data &data, const Event &event, bool recordCounts) {
    if (_reversibleReactions.size() != 1) {
        throw std::runtime_error("Currently works only for one reversible reaction");
    }
    const auto &ctx = kernel->context();
    const auto box = kernel->context().boxSize().data();
    const auto pbc = kernel->context().periodicBoundaryConditions().data();
    auto &model = kernel->getSCPUKernelStateModel();
    scpu_data::new_entries newParticles{};
    std::vector<scpu_data::entry_index> decayedEntries{};

    scalar energyDelta = 0;
    if (event.nEducts == 1) {
        auto &entry1 = data.entry_at(event.idx1);
        auto reaction = ctx.reactions().order1ByType(event.t1)[event.reactionIndex];
        if (reaction->type() != readdy::model::reactions::ReactionType::Fission) {
            throw std::runtime_error("Detailed balance reaction handler only handles Fission or Fusion");
        }
        auto n3 = readdy::model::rnd::normal3<readdy::scalar>(0, 1);
        n3 /= std::sqrt(n3 * n3);
        auto revReaction = _reversibleReactions.front(); // @fixme search for reaction instead of using the first
        auto distance = revReaction.drawFissionDistance();
        Vec3 difference(distance, 0, 0); // orientation does not matter for energy
        scalar energyGain = 0.;
        for (const auto &p : revReaction.lhsPotentials) {
            energyGain += p->calculateEnergy(difference);
        }
        energyDelta = energyGain;
        // create new particles, do not re-use entries, such that all 'old' particles end up in the decayedEentries,
        // which is required for constructing the appropriate backward update
        readdy::model::Particle p1(entry1.position() - reaction->weight2() * distance * n3,
                                  reaction->products()[1]);
        bcs::fixPosition(p1.getPos(), box, pbc);
        newParticles.emplace_back(p1);

        readdy::model::Particle p2(entry1.position() + reaction->weight1() * distance * n3,
                                   reaction->products()[0]);
        bcs::fixPosition(p2.getPos(), box, pbc);
        newParticles.emplace_back(p2);

        decayedEntries.push_back(event.idx1);
        if (recordCounts) {
            model.reactionCounts().at(reaction->id())++;
        }
    } else if (event.nEducts == 2) {
        auto& entry1 = data.entry_at(event.idx1);
        auto& entry2 = data.entry_at(event.idx2);
        auto reaction = ctx.reactions().order2ByType(event.t1, event.t2)[event.reactionIndex];
        if (reaction->type() != readdy::model::reactions::ReactionType::Fusion) {
            throw std::runtime_error("Detailed balance reaction handler only handles Fission or Fusion");
        }
        const auto e1Pos = entry1.pos;
        const auto e2Pos = entry2.pos;
        auto revReaction = _reversibleReactions.front();
        const auto difference = bcs::shortestDifference(e1Pos, e2Pos, box, pbc);
        scalar energyLoss = 0.;
        for (const auto &p : revReaction.lhsPotentials) {
            energyLoss += p->calculateEnergy(difference);
        }
        energyDelta = -1. * energyLoss;
        Vec3 position;
        if (reaction->educts()[0] == entry1.type) {
            position = entry1.pos + reaction->weight1() * difference;
        } else {
            position = entry1.pos + reaction->weight2() * difference;
        }
        readdy::model::Particle particle(position, reaction->products()[0]);
        bcs::fixPosition(particle.getPos(), box, pbc);
        newParticles.emplace_back(particle);
        decayedEntries.push_back(event.idx1);
        decayedEntries.push_back(event.idx2);
        if (recordCounts) {
            model.reactionCounts().at(reaction->id())++;
        }
    } else {
        throw std::runtime_error("This should not happen");
    }
    log::trace("DB performEvent energyDelta {}, newParticles.size() {}, decayedEntries.size() {}",
               energyDelta, newParticles.size(), decayedEntries.size());
    return std::make_pair(std::make_pair(std::move(newParticles), decayedEntries), energyDelta);
}

void SCPUDetailedBalance::calculateForcesEnergies() {
    const auto &context = kernel->context();

    auto &stateModel = kernel->getSCPUKernelStateModel();
    auto &data = *stateModel.getParticleData();
    auto &neighborList = *stateModel.getNeighborList();

    stateModel.energy() = 0;
    stateModel.virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};

    const auto &potentials = context.potentials();
    auto &topologies = stateModel.topologies();
    if (!potentials.potentialsOrder1().empty() || !potentials.potentialsOrder2().empty() || !topologies.empty()) {
        std::for_each(data.begin(), data.end(), [](auto &entry) {
            entry.force = {0, 0, 0};
        });
    }

    // update forces and energy order 1 potentials
    if (!potentials.potentialsOrder1().empty()) {
        {
            std::transform(data.begin(), data.end(), data.begin(), [&potentials, &stateModel](auto &entry) {
                if (!entry.deactivated) {
                    for (const auto &po1 : potentials.potentialsOf(entry.type)) {
                        po1->calculateForceAndEnergy(entry.force, stateModel.energy(), entry.position());
                    }
                }
                return entry;
            });
        }
    }

    // update forces and energy order 2 potentials
    if (!potentials.potentialsOrder2().empty()) {
        const auto &box = context.boxSize().data();
        const auto &pbc = context.periodicBoundaryConditions().data();

        for (auto cell = 0_z; cell < neighborList.nCells(); ++cell) {
            for (auto it = neighborList.particlesBegin(cell); it != neighborList.particlesEnd(cell); ++it) {
                auto pidx = *it;
                auto &entry = data.entry_at(pidx);
                const auto &pots = potentials.potentialsOrder2(entry.type);
                neighborList.forEachNeighbor(it, cell, [&](const std::size_t neighbor) {
                    auto &neighborEntry = data.entry_at(neighbor);
                    auto itPot = pots.find(neighborEntry.type);
                    if (itPot != pots.end()) {
                        Vec3 forceVec{0, 0, 0};
                        auto x_ij = bcs::shortestDifference(entry.position(), neighborEntry.position(), box, pbc);
                        for (const auto &potential : itPot->second) {
                            potential->calculateForceAndEnergy(forceVec, stateModel.energy(), x_ij);
                        }
                        entry.force += forceVec;
                        neighborEntry.force -= forceVec;
                    }
                });
            }
        }
    }
}

}
}
}
}
}


