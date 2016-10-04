/**
 * << detailed description >>
 *
 * @file SingleCPUDefaultReactionProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>

using _rdy_particle_t = readdy::model::Particle;

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace programs {

namespace reactions {
UncontrolledApproximation::UncontrolledApproximation(SingleCPUKernel const *const kernel)
        : kernel(kernel) {
}

void UncontrolledApproximation::execute() {
    const auto &ctx = kernel->getKernelContext();
    const auto &dist = ctx.getDistSquaredFun();
    const auto &fixPos = ctx.getFixPositionFun();
    const auto &dt = ctx.getTimeStep();
    auto data = kernel->getKernelStateModel().getParticleData();
    auto rnd = readdy::model::RandomProvider();
    std::vector<_rdy_particle_t> newParticles{};
    std::vector<std::function<void()>> events{};

    // reactions with one educt
    {
        auto it_type = data->begin_types();

        while (it_type != data->end_types()) {
            // gather reactions
            const auto &reactions = ctx.getOrder1Reactions(*it_type);
            for (const auto &reaction : reactions) {
                auto r = reaction->getRate() * dt;
                if (rnd.getUniform() < r) {
                    const size_t particleIdx = (const size_t) (it_type - data->begin_types());
                    events.push_back([particleIdx, &newParticles, &reaction, this] {
                        auto &&_data = kernel->getKernelStateModel().getParticleData();
                        if (_data->isMarkedForDeactivation(particleIdx)) return;
                        _data->markForDeactivation(particleIdx);

                        switch (reaction->getNProducts()) {
                            case 0: {
                                // no operation, just deactivation
                                // (read out loud with saxony accent)
                                break;
                            }
                            case 1: {
                                if (mapping_11.find(reaction->getId()) != mapping_11.end()) {
                                    newParticles.push_back(
                                            mapping_11[reaction->getId()]((*_data)[particleIdx])
                                    );
                                } else {
                                    const auto particle = (*_data)[particleIdx];
                                    _rdy_particle_t outParticle1{};
                                    reaction->perform(particle, particle, outParticle1,
                                                      outParticle1);
                                    newParticles.push_back(outParticle1);
                                }
                                break;
                            }
                            case 2: {
                                const auto particle = (*_data)[particleIdx];
                                _rdy_particle_t outParticle1{}, outParticle2{};
                                if (mapping_12.find(reaction->getId()) != mapping_12.end()) {
                                    mapping_12[reaction->getId()](particle, outParticle1,
                                                                  outParticle2);
                                } else {
                                    reaction->perform(particle, particle, outParticle1,
                                                      outParticle2);
                                }
                                newParticles.push_back(outParticle1);
                                newParticles.push_back(outParticle2);
                                break;
                            }
                            default: {
                                BOOST_LOG_TRIVIAL(error) << "This should not happen!";
                            }
                        }
                    });
                }
            }
            ++it_type;
        }
    }

    // reactions with two educts
    {
        const auto *neighborList = kernel->getKernelStateModel().getNeighborList();
        for (auto &&it = neighborList->begin(); it != neighborList->end(); ++it) {
            const auto idx1 = it->idx1, idx2 = it->idx2;
            const auto &reactions = ctx.getOrder2Reactions(
                    *(data->begin_types() + idx1), *(data->begin_types() + idx2)
            );

            const auto distSquared = dist(
                    *(data->begin_positions() + idx1), *(data->begin_positions() + idx2)
            );

            for (const auto &reaction : reactions) {
                // if close enough and coin flip successful
                if (distSquared < reaction->getEductDistance() * reaction->getEductDistance()
                    && rnd.getUniform() < reaction->getRate() * dt) {
                    events.push_back([idx1, idx2, this, &newParticles, &reaction] {
                        auto &&_data = kernel->getKernelStateModel().getParticleData();
                        if (_data->isMarkedForDeactivation(idx1)) return;
                        if (_data->isMarkedForDeactivation(idx2)) return;
                        _data->markForDeactivation(idx1);
                        _data->markForDeactivation(idx2);
                        const auto inParticle1 = (*_data)[idx1];
                        const auto inParticle2 = (*_data)[idx2];
                        switch (reaction->getNProducts()) {
                            case 1: {
                                if (mapping_21.find(reaction->getId()) != mapping_21.end()) {
                                    newParticles.push_back(
                                            mapping_21[reaction->getId()](inParticle1,
                                                                          inParticle2));
                                } else {
                                    _rdy_particle_t outParticle1{};
                                    reaction->perform(inParticle1, inParticle2, outParticle1,
                                                      outParticle1);
                                    newParticles.push_back(outParticle1);
                                }
                                break;
                            }
                            case 2: {
                                _rdy_particle_t outParticle1{}, outParticle2{};
                                if (mapping_22.find(reaction->getId()) != mapping_22.end()) {
                                    mapping_22[reaction->getId()](inParticle1, inParticle2,
                                                                  outParticle1, outParticle2);
                                } else {
                                    reaction->perform(inParticle1, inParticle2, outParticle1,
                                                      outParticle2);
                                }
                                newParticles.push_back(outParticle1);
                                newParticles.push_back(outParticle2);
                                break;
                            }
                            default: {
                                BOOST_LOG_TRIVIAL(error) << "This should not happen!";
                            }
                        }
                    });
                }
            }
        }
    }
    // shuffle reactions
    std::random_shuffle(events.begin(), events.end());

    // execute reactions
    std::for_each(events.begin(), events.end(), [](const std::function<void()> &f) { f(); });

    // reposition particles to respect the periodic b.c.
    std::for_each(newParticles.begin(), newParticles.end(),
                  [&fixPos](_rdy_particle_t &p) { fixPos(p.getPos()); });

    // update data structure
    data->deactivateMarked();
    data->addParticles(newParticles);
}

void UncontrolledApproximation::registerReactionScheme_11(const std::string &reactionName,
                                                          reaction_11 fun) {
    auto reaction = kernel->getKernelContext().getReactionOrder1WithName(reactionName);
    if (reaction) {
        if (reaction->getNEducts() == 1 && reaction->getNProducts() == 1) {
            mapping_11.emplace(reaction->getId(), fun);
        } else {
            throw std::runtime_error(
                    "Reaction did not have exactly one product and exactly one educt.");
        }
    } else {
        throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
    }
}

void UncontrolledApproximation::registerReactionScheme_12(const std::string &reactionName,
                                                          reaction_12 fun) {
    auto reaction = kernel->getKernelContext().getReactionOrder1WithName(reactionName);
    if (reaction) {
        if (reaction->getNEducts() == 1 && reaction->getNProducts() == 2) {
            mapping_12.emplace(reaction->getId(), fun);
        } else {
            throw std::runtime_error(
                    "Reaction did not have exactly two products and exactly one educt.");
        }
    } else {
        throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
    }
}

void UncontrolledApproximation::registerReactionScheme_21(const std::string &reactionName,
                                                          reaction_21 fun) {
    auto reaction = kernel->getKernelContext().getReactionOrder2WithName(reactionName);
    if (reaction) {
        if (reaction->getNEducts() == 2 && reaction->getNProducts() == 1) {
            mapping_21.emplace(reaction->getId(), fun);
        } else {
            throw std::runtime_error(
                    "Reaction did not have exactly one product and exactly two educts.");
        }
    } else {
        throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
    }
}

void UncontrolledApproximation::registerReactionScheme_22(const std::string &reactionName,
                                                          reaction_22 fun) {
    auto reaction = kernel->getKernelContext().getReactionOrder2WithName(reactionName);
    if (reaction) {
        if (reaction->getNEducts() == 2 && reaction->getNProducts() == 2) {
            mapping_22.emplace(reaction->getId(), fun);
        } else {
            throw std::runtime_error(
                    "Reaction did not have exactly two products and exactly two educts.");
        }
    } else {
        throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
    }
}


ReactionEvent::ReactionEvent(unsigned int nEducts, unsigned int nProducts,
                             index_type idx1, index_type idx2,
                             double reactionRate, double cumulativeRate,
                             index_type reactionIdx, unsigned int t1, unsigned int t2)
        : nEducts(nEducts), nProducts(nProducts), idx1(idx1), idx2(idx2), reactionRate(reactionRate),
          cumulativeRate(cumulativeRate), reactionIdx(reactionIdx), t1(t1), t2(t2) {}

std::ostream &operator<<(std::ostream &os, const ReactionEvent &evt) {
    os << "ReactionEvent(" << evt.idx1 << "[type=" << evt.t1 << "]";
    if (evt.nEducts == 2) {
        os << " + " << evt.idx2 << "[type=" << evt.t2 << "]";
    }
    os << ", rate=" << evt.reactionRate << ", cumulativeRate=" << evt.cumulativeRate
       << ", reactionIdx=" << evt.reactionIdx;
    return os;
}


std::vector<readdy::model::Particle> Gillespie::handleEvents(std::vector<ReactionEvent> events, double alpha) {
    using _rdy_particle_t = readdy::model::Particle;
    std::vector<_rdy_particle_t> newParticles{};

    const auto &ctx = kernel->getKernelContext();
    auto rnd = readdy::model::RandomProvider();
    auto data = kernel->getKernelStateModel().getParticleData();
    const auto &dt = ctx.getTimeStep();
    /**
     * Handle gathered reaction events
     */
    {
        std::size_t nDeactivated = 0;
        const std::size_t nEvents = events.size();
        while (nDeactivated < nEvents) {
            alpha = (*(events.end() - nDeactivated - 1)).cumulativeRate;
            const auto x = rnd.getUniform(0, alpha);
            const auto eventIt = std::lower_bound(
                    events.begin(), events.end() - nDeactivated, x,
                    [](const ReactionEvent &elem1, double elem2) {
                        return elem1.cumulativeRate < elem2;
                    }
            );
            const auto event = *eventIt;
            if (eventIt == events.end() - nDeactivated) {
                throw std::runtime_error("this should not happen (event not found)");
            }
            if (rnd.getUniform() < event.reactionRate * dt) {
                /**
                 * Perform reaction
                 */
                {
                    const auto p1 = data->operator[](event.idx1);
                    _rdy_particle_t pOut1{}, pOut2{};
                    if (event.nEducts == 1) {
                        auto reaction = ctx.getOrder1Reactions(event.t1)[event.reactionIdx];
                        if (reaction->getNProducts() == 1) {
                            reaction->perform(p1, p1, pOut1, pOut2);
                            newParticles.push_back(pOut1);
                        } else if (reaction->getNProducts() == 2) {
                            reaction->perform(p1, data->operator[](event.idx2), pOut1, pOut2);
                            newParticles.push_back(pOut1);
                            newParticles.push_back(pOut2);
                        }
                    } else {
                        auto reaction = ctx.getOrder2Reactions(event.t1,
                                                               event.t2)[event.reactionIdx];
                        const auto p2 = data->operator[](event.idx2);
                        if (reaction->getNProducts() == 1) {
                            reaction->perform(p1, p2, pOut1, pOut2);
                            newParticles.push_back(pOut1);
                        } else if (reaction->getNProducts() == 2) {
                            reaction->perform(p1, p2, pOut1, pOut2);
                            newParticles.push_back(pOut1);
                            newParticles.push_back(pOut2);
                        }
                    }
                }
                /**
                 * deactivate events whose educts have disappeared (including the just handled one)
                 */
                {
                    auto _it = events.begin();
                    double cumsum = 0.0;
                    const auto idx1 = event.idx1;
                    if (event.nEducts == 1) {
                        data->markForDeactivation((size_t) idx1);
                        while (_it < events.end() - nDeactivated) {
                            if ((*_it).idx1 == idx1 ||
                                ((*_it).nEducts == 2 && (*_it).idx2 == idx1)) {
                                nDeactivated++;
                                std::iter_swap(_it, events.end() - nDeactivated);
                                cumsum += (*_it).reactionRate;
                                (*_it).cumulativeRate = cumsum;
                            } else {
                                ++_it;
                            }
                        }
                    } else {
                        const auto idx2 = event.idx2;
                        data->markForDeactivation((size_t) event.idx1);
                        data->markForDeactivation((size_t) event.idx2);
                        while (_it < events.end() - nDeactivated) {
                            if ((*_it).idx1 == idx1 || (*_it).idx1 == idx2
                                ||
                                ((*_it).nEducts == 2 &&
                                 ((*_it).idx2 == idx1 || (*_it).idx2 == idx2))) {
                                nDeactivated++;
                                std::iter_swap(_it, events.end() - nDeactivated);
                                cumsum += (*_it).reactionRate;
                                (*_it).cumulativeRate = cumsum;
                            } else {
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
    return newParticles;
}

std::vector<ReactionEvent> Gillespie::gatherEvents(double &alpha) {
    std::vector<ReactionEvent> events;
    const auto &ctx = kernel->getKernelContext();
    auto data = kernel->getKernelStateModel().getParticleData();
    const auto &dist = ctx.getDistSquaredFun();
    /**
    * Reactions with one educt
    */
    {
        auto it_type = data->begin_types();
        const auto end = data->end_types();
        while (it_type != end) {
            const auto &reactions = ctx.getOrder1Reactions(*it_type);
            for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                const auto rate = (*it)->getRate();
                if (rate > 0) {
                    alpha += rate;
                    events.push_back(
                            {1, (*it)->getNProducts(), (_reaction_idx_t) (it_type - data->begin_types()), 0, rate,
                             alpha,
                             (_reaction_idx_t) (it - reactions.begin()), *it_type, 0});
                }
            }
            ++it_type;
        }
    }

    /**
     * Reactions with two educts
     */
    {
        const auto *neighborList = kernel->getKernelStateModel().getNeighborList();
        auto typesBegin = data->begin_types();
        for (auto &&nl_it = neighborList->begin(); nl_it != neighborList->end(); ++nl_it) {
            const _reaction_idx_t idx1 = nl_it->idx1, idx2 = nl_it->idx2;
            const auto &reactions = ctx.getOrder2Reactions(
                    *(data->begin_types() + idx1), *(data->begin_types() + idx2)
            );

            const auto distSquared = dist(
                    *(data->begin_positions() + idx1), *(data->begin_positions() + idx2)
            );

            for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                // if close enough
                const auto reaction = *it;
                if (distSquared < reaction->getEductDistance() * reaction->getEductDistance()) {
                    const auto rate = reaction->getRate();
                    if (rate > 0) {
                        alpha += rate;
                        events.push_back(
                                {2, reaction->getNProducts(), idx1, idx2, rate, alpha,
                                 (_reaction_idx_t) (it - reactions.begin()),
                                 *(typesBegin + idx1), *(typesBegin + idx2)});
                    }
                }
            }
        }
    }
    return events;
}
}
}
}
}
}


