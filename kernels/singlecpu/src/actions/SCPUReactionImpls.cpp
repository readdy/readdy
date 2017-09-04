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
          cumulativeRate(cumulativeRate), reactionIdx(reactionIdx), t1(t1), t2(t2) {
}

std::ostream &operator<<(std::ostream &os, const Event &evt) {
    os << "Event(" << evt.idx1 << "[type=" << evt.t1 << "]";
    if (evt.nEducts == 2) {
        os << " + " << evt.idx2 << "[type=" << evt.t2 << "]";
    }
    os << ", rate=" << evt.reactionRate << ", cumulativeRate=" << evt.cumulativeRate
       << ", reactionIdx=" << evt.reactionIdx;
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

std::vector<event_t> findEvents(const SCPUKernel *const kernel, scalar dt, bool approximateRate = true) {
    std::vector<event_t> eventsUpdate;
    auto &stateModel = kernel->getSCPUKernelStateModel();
    auto &data = *stateModel.getParticleData();
    const auto &d2 = kernel->getKernelContext().getDistSquaredFun();
    {
        std::size_t idx = 0;
        for (const auto &e : data) {
            if (!e.is_deactivated()) {
                // order 1
                const auto &reactions = kernel->getKernelContext().reactions().order1_by_type(e.type);
                for (auto it_reactions = reactions.begin(); it_reactions != reactions.end(); ++it_reactions) {
                    const auto rate = (*it_reactions)->getRate();
                    if (rate > 0 && shouldPerformEvent(rate, dt, approximateRate)) {
                        //unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, scalar reactionRate,
                        //scalar cumulativeRate, reaction_index_type reactionIdx, particletype_t t1, particletype_t t2
                        Event evt{1, (*it_reactions)->getNProducts(), idx, idx, rate, 0,
                                  static_cast<std::size_t>(it_reactions - reactions.begin()),
                                  e.type, 0};
                        eventsUpdate.push_back(evt);
                    }
                }
            }
            ++idx;
        }
    }
    auto it_nl = stateModel.getNeighborList()->cbegin();
    for (; it_nl != stateModel.getNeighborList()->cend(); ++it_nl) {
        const auto &entry = data.entry_at(it_nl->idx1);
        if (!entry.is_deactivated()) {

            // order 2
            const auto &neighbor = data.entry_at(it_nl->idx2);
            const auto &reactions = kernel->getKernelContext().reactions().order2_by_type(entry.type, neighbor.type);
            if (!reactions.empty()) {
                const auto distSquared = d2(neighbor.position(), entry.position());
                for (auto it_reactions = reactions.begin(); it_reactions < reactions.end(); ++it_reactions) {
                    const auto &react = *it_reactions;
                    const auto rate = react->getRate();
                    if (rate > 0 && distSquared < react->getEductDistanceSquared()
                        && shouldPerformEvent(rate, dt, approximateRate)) {
                        const auto reaction_index = static_cast<event_t::reaction_index_type>(it_reactions -
                                                                                              reactions.begin());
                        eventsUpdate.emplace_back(2, react->getNProducts(), it_nl->idx1, it_nl->idx2, rate, 0, reaction_index,
                                 entry.type, neighbor.type);
                    }
                }
            }
        }
    }
    return eventsUpdate;
}

void SCPUUncontrolledApproximation::perform(util::PerformanceNode &node) {
    auto timer = node.timeit();
    const auto &ctx = kernel->getKernelContext();
    const auto &fixPos = ctx.getFixPositionFun();
    SCPUStateModel &stateModel = kernel->getSCPUKernelStateModel();
    if(ctx.recordReactionsWithPositions()) stateModel.reactionRecords().clear();
    if(ctx.recordReactionCounts()) {
        readdy::model::observables::ReactionCounts::initializeCounts(stateModel.reactionCounts(), ctx);
    }
    auto &data = *stateModel.getParticleData();
    auto events = findEvents(kernel, timeStep, true);

    // shuffle reactions
    std::shuffle(events.begin(), events.end(), std::mt19937(std::random_device()()));

    // execute reactions
    {
        scpu_data::new_entries newParticles{};
        std::vector<scpu_data::entry_index> decayedEntries{};

        for (auto it = events.begin(); it != events.end(); ++it) {
            auto &event = *it;
            if (event.cumulativeRate == 0) {
                auto entry1 = event.idx1;
                if (event.nEducts == 1) {
                    auto reaction = ctx.reactions().order1_by_type(event.t1)[event.reactionIdx];
                    if(ctx.recordReactionsWithPositions()) {
                        reaction_record record;
                        record.reactionIndex = event.reactionIdx;
                        performReaction(data, entry1, entry1, newParticles, decayedEntries, reaction, fixPos, &record);
                        stateModel.reactionRecords().push_back(record);
                    } else {
                        performReaction(data, entry1, entry1, newParticles, decayedEntries, reaction, fixPos, nullptr);
                    }
                    if(ctx.recordReactionCounts()) {
                        auto &countsOrder1 = std::get<0>(stateModel.reactionCounts());
                        countsOrder1.at(event.t1).at(event.reactionIdx)++;
                    }
                    for (auto _it2 = it + 1; _it2 != events.end(); ++_it2) {
                        if (_it2->idx1 == entry1 || _it2->idx2 == entry1) {
                            _it2->cumulativeRate = 1;
                        }
                    }
                } else {
                    auto reaction = ctx.reactions().order2_by_type(event.t1, event.t2)[event.reactionIdx];
                    if(ctx.recordReactionsWithPositions()) {
                        reaction_record record;
                        record.reactionIndex = event.reactionIdx;
                        performReaction(data, entry1, event.idx2, newParticles, decayedEntries, reaction, fixPos, &record);
                        stateModel.reactionRecords().push_back(record);
                    } else {
                        performReaction(data, entry1, event.idx2, newParticles, decayedEntries, reaction, fixPos, nullptr);
                    }
                    if(ctx.recordReactionCounts()) {
                        auto &countsOrder2 = std::get<1>(stateModel.reactionCounts());
                        countsOrder2.at(std::tie(event.t1, event.t2)).at(event.reactionIdx)++;
                    }
                    for (auto _it2 = it + 1; _it2 != events.end(); ++_it2) {
                        if (_it2->idx1 == entry1 || _it2->idx2 == entry1 ||
                            _it2->idx1 == event.idx2 || _it2->idx2 == event.idx2) {
                            _it2->cumulativeRate = 1;
                        }
                    }
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
                  const readdy::kernel::scpu::model::SCPUNeighborList &nl, const scpu_data &data, scalar &alpha,
                  std::vector<event_t> &events, const readdy::model::KernelContext::dist_squared_fun &d2) {
    {
        std::size_t index = 0;
        for (const auto &entry : data) {
            if(!entry.is_deactivated()) {
                const auto &reactions = kernel->getKernelContext().reactions().order1_by_type(entry.type);
                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    const auto rate = (*it)->getRate();
                    if (rate > 0) {
                        alpha += rate;
                        events.emplace_back(1, (*it)->getNProducts(), index, 0, rate, alpha,
                                 static_cast<event_t::reaction_index_type>(it - reactions.begin()),
                                 entry.type, 0);
                    }
                }
            }
            ++index;
        }
    }
    for (auto it_nl : nl) {
        auto &entry = data.entry_at(it_nl.idx1);
        if (!entry.is_deactivated()) {
            // order 2
            const auto &neighbor = data.entry_at(it_nl.idx2);
            const auto &reactions = kernel->getKernelContext().reactions().order2_by_type(entry.type,
                                                                                                 neighbor.type);
            if (!reactions.empty()) {
                const auto distSquared = d2(neighbor.position(), entry.position());
                for (auto it = reactions.begin(); it < reactions.end(); ++it) {
                    const auto &react = *it;
                    const auto rate = react->getRate();
                    if (rate > 0 && distSquared < react->getEductDistanceSquared()) {
                        alpha += rate;
                        events.emplace_back(2, react->getNProducts(), it_nl.idx1, it_nl.idx2,
                                          rate, alpha,
                                          static_cast<event_t::reaction_index_type>(it - reactions.begin()),
                                          entry.type, neighbor.type);
                    }
                }
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
        const auto &ctx = kernel->getKernelContext();
        const auto &fixPos = ctx.getFixPositionFun();
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
                            auto reaction = ctx.reactions().order1_by_type(event.t1)[event.reactionIdx];
                            if(ctx.recordReactionsWithPositions()) {
                                reaction_record record;
                                record.reactionIndex = event.reactionIdx;
                                performReaction(*data, entry1, entry1, newParticles, decayedEntries, reaction, fixPos,
                                                &record);
                                model.reactionRecords().push_back(std::move(record));
                            } else {
                                performReaction(*data, entry1, entry1, newParticles, decayedEntries, reaction, fixPos,
                                                nullptr);
                            }
                            if(ctx.recordReactionCounts()) {
                                auto &countsOrder1 = std::get<0>(model.reactionCounts());
                                countsOrder1.at(event.t1).at(event.reactionIdx)++;
                            }
                        } else {
                            auto reaction = ctx.reactions().order2_by_type(event.t1, event.t2)[event.reactionIdx];
                            if(ctx.recordReactionsWithPositions()) {
                                reaction_record record;
                                record.reactionIndex = event.reactionIdx;
                                performReaction(*data, entry1, event.idx2, newParticles, decayedEntries, reaction, fixPos, &record);
                                model.reactionRecords().push_back(std::move(record));
                            } else {
                                performReaction(*data, entry1, event.idx2, newParticles, decayedEntries, reaction, fixPos, nullptr);
                            }
                            if(ctx.recordReactionCounts()) {
                                auto &countsOrder2 = std::get<1>(model.reactionCounts());
                                countsOrder2.at(std::tie(event.t1, event.t2)).at(event.reactionIdx)++;
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

void SCPUGillespie::perform(util::PerformanceNode &node) {
    auto timer = node.timeit();
    const auto &ctx = kernel->getKernelContext();
    if(ctx.reactions().n_order1() == 0 && ctx.reactions().n_order2() == 0) {
        return;
    }
    auto &stateModel = kernel->getSCPUKernelStateModel();
    if(ctx.recordReactionsWithPositions()) stateModel.reactionRecords().clear();
    if(ctx.recordReactionCounts()) {
        readdy::model::observables::ReactionCounts::initializeCounts(stateModel.reactionCounts(), ctx);
    }
    auto data = stateModel.getParticleData();
    const auto &dist = ctx.getDistSquaredFun();
    const auto &fixPos = ctx.getFixPositionFun();
    const auto nl = stateModel.getNeighborList();

    scalar alpha = 0.0;
    std::vector<event_t> events;
    gatherEvents(kernel, *nl, *data, alpha, events, dist);
    auto particlesUpdate = handleEventsGillespie(kernel, timeStep, false, true, std::move(events));

    // update data structure
    data->update(std::move(particlesUpdate));
}
}
}
}
}
}


