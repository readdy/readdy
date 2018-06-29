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
 * @file SCPUReactionImpls.cpp
 * @brief Implementations of reaction algorithms for the SingleCPU kernel
 * @author clonker
 * @author chrisfroe
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
             scalar cumulativeRate, reaction_index_type reactionIdx, ParticleTypeId t1, ParticleTypeId t2)
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
    auto events = findEvents(kernel, timeStep(), false);

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
    auto particlesUpdate = handleEventsGillespie(kernel, timeStep(), false, false, std::move(events));

    // update data structure
    data->update(std::move(particlesUpdate));
}

ParticleBackup::ParticleBackup(event_t event,
                               const readdy::model::actions::reactions::ReversibleReactionConfig *revReaction,
                               const readdy::model::reactions::Reaction *reaction, const scpu_data *data) {
    switch (revReaction->reversibleType) {
        case readdy::model::actions::reactions::FusionFission: {
            nParticles = event.nEducts;
            idx1 = event.idx1;
            t1 = event.t1;
            idx2 = event.idx2;
            t2 = event.t2;
            pos1 = data->entry_at(idx1).position();
            pos2 = data->entry_at(idx2).position();
            break;
        }
        case readdy::model::actions::reactions::ConversionConversion: {
            nParticles = 1;
            idx1 = event.idx1;
            t1 = event.t1;
            pos1 = data->entry_at(idx1).position();
            break;
        }
        case readdy::model::actions::reactions::EnzymaticEnzymatic:{
            nParticles = 1;
            // find out which particle is the catalyst in A + C -> B + C
            if (event.t1 == reaction->educts()[1]) {
                idx1 = event.idx2;
                t1 = event.t2;
                pos1 = data->entry_at(idx1).position();
            } else if (event.t2 == reaction->educts()[1]){
                idx1 = event.idx1;
                t1 = event.t1;
                pos1 = data->entry_at(idx1).position();
            } else {
                throw std::runtime_error(
                        fmt::format("None of the event's particles could be identified as catalyst, method: {} file: {}",
                                    "ParticleBackup::ParticleBackup", "SCPUReactionImpls.cpp"));
            }
            break;
        }
        default:
            throw std::runtime_error(fmt::format("Unknown type of reversible reaction, method: {} file: {}",
                                                 "ParticleBackup::ParticleBackup", "SCPUReactionImpls.cpp"));
    }
}

void SCPUDetailedBalance::perform(const readdy::util::PerformanceNode &node) {
    auto t = node.timeit();

    const auto &ctx = kernel->context();
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

    const auto approximateRate = false;

    {
        auto shouldEval = [&](const event_t &event) {
            return shouldPerformEvent(event.rate, timeStep(), approximateRate);
        };

        auto depending = [&](const event_t &e1, const event_t &e2) {
            return (e1.idx1 == e2.idx1 || (e2.nEducts == 2 && e1.idx1 == e2.idx2)
                    || (e1.nEducts == 2 && (e1.idx2 == e2.idx1 || (e2.nEducts == 2 && e1.idx2 == e2.idx2))));
        };

        auto eval = [&](const event_t &event) {
            const readdy::model::reactions::Reaction *reaction;
            const readdy::model::actions::reactions::ReversibleReactionConfig *revReaction;
            std::tie(revReaction, reaction) = findReversibleReaction(event);
            bool isReversibleReaction = (revReaction!=nullptr);

            if (isReversibleReaction) {
                // Perform detailed balance method for reversible reaction event

                const auto particleBackup = ParticleBackup(event, revReaction, reaction, data);
                const auto energyBefore = stateModel.energy();

                reaction_record record;
                model::SCPUParticleData::entries_update forwardUpdate;
                scalar interactionEnergy; // only relevant for FusionFission
                if (ctx.recordReactionsWithPositions()) {
                    std::tie(forwardUpdate, interactionEnergy) = performReversibleReactionEvent(event, revReaction,
                                                                                                reaction, &record);
                } else {
                    std::tie(forwardUpdate, interactionEnergy) = performReversibleReactionEvent(event, revReaction,
                                                                                                reaction, nullptr);
                }
                const auto updateRecord = data->update(std::move(forwardUpdate));
                auto backwardUpdate = generateBackwardUpdate(particleBackup, updateRecord);
                stateModel.updateNeighborList();
                calculateEnergies(node);

                scalar boltzmannFactor = 1.;
                scalar prefactor = 1.;
                switch (revReaction->reversibleType) {
                    case readdy::model::actions::reactions::FusionFission: {
                        // the sign of interactionEnergy was determined automatically in perform...()
                        // thus not checking if forward or backward here
                        prefactor = 1.;
                        boltzmannFactor = std::exp(
                                -1. / ctx.kBT() * ((stateModel.energy() - interactionEnergy) - energyBefore));
                        break;
                    }
                    case readdy::model::actions::reactions::ConversionConversion: {
                        prefactor = 1.;
                        boltzmannFactor = std::exp(-1. / ctx.kBT() * (stateModel.energy() - energyBefore));
                        break;
                    }
                    case readdy::model::actions::reactions::EnzymaticEnzymatic: {
                        if (reaction->educts() == revReaction->lhsTypes) {
                            // forward
                            prefactor = revReaction->acceptancePrefactor;
                        } else if (reaction->educts() == revReaction->lhsTypes) {
                            prefactor = 1. / revReaction->acceptancePrefactor;
                        }
                        boltzmannFactor = std::exp(-1. / ctx.kBT() * (stateModel.energy() - energyBefore));
                        break;
                    }
                    default:
                        throw std::runtime_error(fmt::format("Unknown type of reversible reaction, method: {} file: {}",
                                                             "SCPUDetailedBalance::perform::eval",
                                                             "SCPUReactionImpls.cpp"));
                }
                const scalar acceptance = std::min(1., prefactor * boltzmannFactor);
                log::trace("Acceptance for current event is {}", acceptance);

                if (readdy::model::rnd::uniform_real() < acceptance) {
                    // accept/do-nothing
                    log::trace("accept!");
                    if (ctx.recordReactionsWithPositions()) {
                        stateModel.reactionRecords().push_back(record);
                    }
                    if (ctx.recordReactionCounts()) {
                        stateModel.reactionCounts().at(reaction->id())++;
                    }
                } else {
                    // reject/rollback
                    log::trace("reject! apply backward update");
                    data->update(std::move(backwardUpdate));
                    stateModel.updateNeighborList();
                    calculateEnergies(node);
                    if (energyBefore != stateModel.energy()) {
                        log::warn("reaction move was rejected but energy of state is different "
                                  "after rollback, was {}, is now {}", energyBefore, stateModel.energy());
                    }
                }
            } else {
                // Perform vanilla Doi model with direct update to data structure
                model::SCPUParticleData::entries_update forwardUpdate;

                if (ctx.recordReactionsWithPositions()) {
                    reaction_record record;
                    record.id = reaction->id();
                    performReaction(*data, event.idx1, event.idx2, forwardUpdate.first, forwardUpdate.second, reaction,
                                    ctx, &record);
                    stateModel.reactionRecords().push_back(record);
                } else {
                    performReaction(*data, event.idx1, event.idx2, forwardUpdate.first, forwardUpdate.second, reaction,
                                    ctx, nullptr);
                }

                if(ctx.recordReactionCounts()) {
                    stateModel.reactionCounts().at(reaction->id())++;
                }

                data->update(std::move(forwardUpdate));
                stateModel.updateNeighborList();
                calculateEnergies(node);
            }
        };

        calculateEnergies(node);
        algo::performEvents(events, shouldEval, depending, eval);

    }
}

model::SCPUParticleData::entries_update
SCPUDetailedBalance::generateBackwardUpdate(
        const ParticleBackup &particleBackup,
        const std::vector<model::SCPUParticleData::entry_index> &updateRecord) const {
    // the entries that were created by the forward update
    std::vector<scpu_data::entry_index> decayedEntries = updateRecord;
    scpu_data::new_entries newParticles{};
    if (particleBackup.nParticles == 1) {
        readdy::model::Particle p(particleBackup.pos1, particleBackup.t1);
        newParticles.emplace_back(p);
    } else if (particleBackup.nParticles == 2) {
        readdy::model::Particle p1(particleBackup.pos1, particleBackup.t1);
        readdy::model::Particle p2(particleBackup.pos2, particleBackup.t2);
        newParticles.emplace_back(p1);
        newParticles.emplace_back(p2);
    } else {
        throw std::runtime_error(
                fmt::format("Particle backup can only contain information on one or two particles, method: {} file: {}",
                            "SCPUDetailedBalance::generateBackwardUpdate", "SCPUReactionImpls.cpp"));
    }
    return std::make_pair(std::move(newParticles), decayedEntries);
}

std::pair<model::SCPUParticleData::entries_update, scalar> SCPUDetailedBalance::performReversibleReactionEvent(
        const Event &event, const readdy::model::actions::reactions::ReversibleReactionConfig *reversibleReaction,
        const readdy::model::reactions::Reaction *reaction, reaction_record *record) {
    if (reversibleReaction == nullptr) {
        throw std::runtime_error(fmt::format("Reaction is not reversible, method: {} file: {}",
                                             "SCPUDetailedBalance::performReversibleReactionEvent",
                                             "SCPUReactionImpls.cpp"));
    }

    auto &model = kernel->getSCPUKernelStateModel();
    scpu_data* data = model.getParticleData();
    const auto &ctx = kernel->context();

    const auto box = ctx.boxSize().data();
    const auto pbc = ctx.periodicBoundaryConditions().data();

    scpu_data::new_entries newParticles{};
    std::vector<scpu_data::entry_index> decayedEntries{};

    scalar energyDelta = 0;
    switch (reversibleReaction->reversibleType) {
        case readdy::model::actions::reactions::FusionFission: {
            if (event.nEducts == 1) {
                // backward reaction C --> A + B
                const auto &entry1 = data->entry_at(event.idx1);
                auto n3 = readdy::model::rnd::normal3<readdy::scalar>(0, 1);
                n3 /= std::sqrt(n3 * n3);
                const auto distance = reversibleReaction->drawFissionDistance();
                Vec3 difference(distance, 0, 0); // orientation does not matter for energy
                scalar energyGain = 0.; // calculate U_AB
                for (const auto &p : reversibleReaction->lhsPotentials) {
                    energyGain += p->calculateEnergy(difference);
                }
                energyDelta = energyGain;
                // IMPORTANT
                // create new particles, do not re-use entries,
                // such that all 'old' particles end up in the decayedEentries,
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

                if (record) {
                    record->id = reaction->id();
                    record->type = static_cast<int>(reaction->type());
                    record->where = (entry1.position());
                    bcs::fixPosition(record->where, box, pbc);
                    record->educts[0] = entry1.id;
                    record->educts[1] = entry1.id;
                    record->types_from[0] = entry1.type;
                    record->types_from[1] = entry1.type;
                }
            } else if (event.nEducts == 2) {
                // forward reaction A + B --> C
                const auto &entry1 = data->entry_at(event.idx1);
                const auto &entry2 = data->entry_at(event.idx2);

                const auto e1Pos = entry1.pos;
                const auto e2Pos = entry2.pos;
                const auto difference = bcs::shortestDifference(e1Pos, e2Pos, box, pbc);
                scalar energyLoss = 0.; // calculate U_AB
                for (const auto &p : reversibleReaction->lhsPotentials) {
                    energyLoss += p->calculateEnergy(difference);
                }
                energyDelta = -1. * energyLoss;
                Vec3 position;
                if (reaction->educts()[0] == entry1.type) {
                    position = entry1.pos + reaction->weight1() * difference;
                } else {
                    position = entry1.pos + reaction->weight2() * difference;
                }

                // do not re-use entries for proper construction of backward update
                readdy::model::Particle particle(position, reaction->products()[0]);
                bcs::fixPosition(particle.getPos(), box, pbc);

                newParticles.emplace_back(particle);
                decayedEntries.push_back(event.idx1);
                decayedEntries.push_back(event.idx2);

                if (record) {
                    record->id = reaction->id();
                    record->type = static_cast<int>(reaction->type());
                    record->where = (entry1.position() + entry2.position()) / 2.;
                    bcs::fixPosition(record->where, box, pbc);
                    record->educts[0] = entry1.id;
                    record->educts[1] = entry2.id;
                    record->types_from[0] = entry1.type;
                    record->types_from[1] = entry2.type;
                }
            }
            break;
        }
        case readdy::model::actions::reactions::ConversionConversion: {
            const auto &entry = data->entry_at(event.idx1);

            // do not re-use entries for proper construction of backward update
            readdy::model::Particle particle(entry.position(), reaction->products()[0]);
            newParticles.emplace_back(particle);
            decayedEntries.push_back(event.idx1);
            energyDelta = 0.;

            if(record) {
                record->id = reaction->id();
                record->type = static_cast<int>(reaction->type());
                record->where = (entry.position());
                bcs::fixPosition(record->where, box, pbc);
                record->educts[0] = entry.id;
                record->educts[1] = entry.id;
                record->types_from[0] = entry.type;
                record->types_from[1] = entry.type;
            }

            break;
        }
        case readdy::model::actions::reactions::EnzymaticEnzymatic: {
            // find out which particle is the catalyst in A + C -> B + C
            event_t::index_type catalystIdx;
            event_t::index_type eductIdx;
            if (event.t1 == reaction->educts()[1]) {
                catalystIdx = event.idx1;
                eductIdx = event.idx2;
            } else if (event.t2 == reaction->educts()[1]){
                catalystIdx = event.idx2;
                eductIdx = event.idx1;
            }

            const auto& eductEntry = data->entry_at(eductIdx);
            const auto& catalystEntry = data->entry_at(catalystIdx);

            // do not re-use entries for proper construction of backward update
            readdy::model::Particle particle(eductEntry.position(), reaction->products()[0]);
            newParticles.emplace_back(particle);
            decayedEntries.push_back(eductIdx);
            energyDelta = 0.;

            if(record) {
                record->id = reaction->id();
                record->type = static_cast<int>(reaction->type());
                record->where = (eductEntry.position() + catalystEntry.position()) / 2.;
                bcs::fixPosition(record->where, box, pbc);
                record->educts[0] = eductEntry.id;
                record->educts[1] = catalystEntry.id;
                record->types_from[0] = eductEntry.type;
                record->types_from[1] = catalystEntry.type;
            }

            break;
        }
        default:
            throw std::runtime_error(fmt::format("Unknown type of reversible reaction, method: {} file: {}",
                                                 "SCPUDetailedBalance::performReversibleReactionEvent",
                                                 "SCPUReactionImpls.cpp"));
    }

    return std::make_pair(std::make_pair(std::move(newParticles), decayedEntries), energyDelta);
}

void SCPUDetailedBalance::calculateEnergies(const util::PerformanceNode &node) {
    const auto &context = kernel->context();

    auto &stateModel = kernel->getSCPUKernelStateModel();
    auto &data = *stateModel.getParticleData();
    auto &neighborList = *stateModel.getNeighborList();
    const auto &potentials = context.potentials();

    stateModel.energy() = 0;

    // calculate energies for particles
    auto order1eval = [&](auto &entry){
        for (const auto &po1 : potentials.potentialsOf(entry.type)) {
            stateModel.energy() += po1->calculateEnergy(entry.position());
        }
    };

    // calculate energies for particle pairs
    const auto &box = context.boxSize().data();
    const auto &pbc = context.periodicBoundaryConditions().data();

    auto order2eval = [&](auto &entry, auto &neighborEntry) {
        const auto &pots = potentials.potentialsOrder2(entry.type);
        auto itPot = pots.find(neighborEntry.type);
        if (itPot != std::end(pots)) {
            auto x_ij = bcs::shortestDifference(entry.position(), neighborEntry.position(), box, pbc);
            for (const auto &potential : itPot->second) {
                stateModel.energy() += potential->calculateEnergy(x_ij);
            }
        }
    };

    auto topologyEval = [](auto &topology){ /* noop */ };
    SCPUStateModel::topologies_vec emptyContainer = {};

    algo::evaluateOnContainers(data, order1eval, neighborList, order2eval, emptyContainer, topologyEval, node);
}

std::pair<const readdy::model::actions::reactions::ReversibleReactionConfig *, const readdy::model::reactions::Reaction *>
SCPUDetailedBalance::findReversibleReaction(const Event &event) {
    const auto &ctx = kernel->context();
    if (event.nEducts == 1) {
        const auto &reaction = ctx.reactions().order1ByType(event.t1)[event.reactionIndex];
        auto findIt = _reversibleReactionsMap.find(reaction->id());
        if (findIt != _reversibleReactionsMap.end()) {
            const auto &revReaction = (*findIt).second;
            return std::make_pair(revReaction.get(), reaction);
        } else {
            return std::make_pair(nullptr, reaction);
        }
    } else {
        const auto &reaction = ctx.reactions().order2ByType(event.t1, event.t2)[event.reactionIndex];
        auto findIt = _reversibleReactionsMap.find(reaction->id());
        if (findIt != _reversibleReactionsMap.end()) {
            const auto &revReaction = (*findIt).second;
            return std::make_pair(revReaction.get(), reaction);
        } else {
            return std::make_pair(nullptr, reaction);
        }
    }
}

}
}
}
}
}


