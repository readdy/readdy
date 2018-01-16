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

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {

namespace reactions {

Event::Event(unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, scalar reactionRate,
             scalar cumulativeRate, reaction_index_type reactionIdx, particle_type_type t1, particle_type_type t2)
        : nEducts(nEducts), nProducts(nProducts), idx1(idx1), idx2(idx2), reactionRate(reactionRate),
          cumulativeRate(cumulativeRate), reactionIndex(reactionIdx), t1(t1), t2(t2) {
}

std::ostream &operator<<(std::ostream &os, const Event &evt) {
    os << "Event(" << evt.idx1 << "[type=" << evt.t1 << "]";
    if (evt.nEducts == 2) {
        os << " + " << evt.idx2 << "[type=" << evt.t2 << "]";
    }
    os << ", rate=" << evt.reactionRate << ", cumulativeRate=" << evt.cumulativeRate
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
    const auto &d2 = kernel->context().distSquaredFun();
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
                            const auto distSquared = d2(neighbor.position(), entry.position());
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
                  std::vector<event_t> &events, const readdy::model::Context::dist_squared_fun &d2) {
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
                            const auto distSquared = d2(neighbor.position(), entry.position());
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
            std::size_t nDeactivated = 0;
            const std::size_t nEvents = events.size();
            while (nDeactivated < nEvents) {
                const auto alpha = (*(events.end() - nDeactivated - 1)).cumulativeRate;
                const auto x = readdy::model::rnd::uniform_real(static_cast<scalar>(0.), alpha);
                const auto eventIt = std::lower_bound(
                        events.begin(), events.end() - nDeactivated, x,
                        [](const event_t &elem1, scalar elem2) {
                            return elem1.cumulativeRate < elem2;
                        }
                );
                const auto &event = *eventIt;
                if (eventIt == events.end() - nDeactivated) {
                    throw std::runtime_error("this should not happen (event not found)");
                }
                if (filterEventsInAdvance || shouldPerformEvent(event.reactionRate, timeStep, approximateRate)) {
                    /**
                     * Perform reaction
                     */
                    {
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
                    }
                    /**
                     * deactivate events whose educts have disappeared (including the just handled one)
                     */
                    {
                        auto _it = events.begin();
                        scalar cumsum = 0.0;
                        const auto idx1 = event.idx1;
                        if (event.nEducts == 1) {
                            while (_it < events.end() - nDeactivated) {
                                if (_it->idx1 == idx1 ||
                                    (_it->nEducts == 2 && _it->idx2 == idx1)) {
                                    ++nDeactivated;
                                    std::iter_swap(_it, events.end() - nDeactivated);
                                } else {
                                    cumsum += _it->reactionRate;
                                    _it->cumulativeRate = cumsum;
                                    ++_it;
                                }
                            }
                        } else {
                            const auto idx2 = event.idx2;
                            while (_it < events.end() - nDeactivated) {
                                if (_it->idx1 == idx1 || _it->idx1 == idx2 ||
                                    (_it->nEducts == 2 &&
                                     (_it->idx2 == idx1 || _it->idx2 == idx2))) {
                                    ++nDeactivated;
                                    std::iter_swap(_it, events.end() - nDeactivated);
                                } else {
                                    _it->cumulativeRate = cumsum;
                                    cumsum += _it->reactionRate;
                                    ++_it;
                                }
                            }
                        }

                    }
                } else {
                    nDeactivated++;
                    std::iter_swap(eventIt, events.end() - nDeactivated);
                    (*eventIt).cumulativeRate = (*eventIt).reactionRate;
                    if (eventIt > events.begin()) {
                        (*eventIt).cumulativeRate += (*(eventIt - 1)).cumulativeRate;
                    }
                    auto cumsum = (*eventIt).cumulativeRate;
                    for (auto _it = eventIt + 1; _it < events.end() - nDeactivated; ++_it) {
                        cumsum += (*_it).reactionRate;
                        (*_it).cumulativeRate = cumsum;
                    }
                }
            }
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
    const auto &dist = ctx.distSquaredFun();
    const auto &fixPos = ctx.fixPositionFun();
    const auto nl = stateModel.getNeighborList();

    scalar alpha = 0.0;
    std::vector<event_t> events;
    gatherEvents(kernel, *nl, *data, alpha, events, dist);
    auto particlesUpdate = handleEventsGillespie(kernel, timeStep, false, false, std::move(events));

    // update data structure
    data->update(std::move(particlesUpdate));
}
}
}
}
}
}


