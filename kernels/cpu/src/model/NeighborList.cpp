/**
 * << detailed description >>
 *
 * @file NeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 08.09.16
 */

#include <readdy/kernel/cpu/model/NeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

struct NeighborList::Box {
    std::vector<Box *> neighbors{};
    std::vector<particle_index> particleIndices{};
    const box_index id = 0;
    const bool enoughBoxes;

    // dirty flag indicating whether the box and its neighboring boxes have to be re-created
    bool dirty = true;

    Box(box_index i, box_index j, box_index k, const std::array<box_index, 3> &nBoxes)
            : id(k + j * nBoxes[2] + i * nBoxes[2] * nBoxes[1]),
              enoughBoxes(nBoxes[0] >= 5 && nBoxes[1] >= 5 && nBoxes[2] >= 5) {
    }

    void addNeighbor(Box *box) {
        if (box && box->id != id
            && (enoughBoxes || std::find(neighbors.begin(), neighbors.end(), box) == neighbors.end())) {
            neighbors.push_back(box);
        }
    }

    friend bool operator==(const Box &lhs, const Box &rhs) {
        return lhs.id == rhs.id;
    }

    friend bool operator!=(const Box &lhs, const Box &rhs) {
        return !(lhs == rhs);
    }
};

NeighborList::NeighborList(const readdy::model::KernelContext *const context, util::Config const *const config)
        : ctx(context), config(config), boxes({}), simBoxSize(ctx->getBoxSize()) {}

void NeighborList::setupBoxes() {
    if (boxes.empty()) {
        double maxCutoff = 0;
        for (auto &&e : ctx->getAllOrder2RegisteredPotentialTypes()) {
            for (auto &&p : ctx->getOrder2Potentials(std::get<0>(e), std::get<1>(e))) {
                maxCutoff = maxCutoff < p->getCutoffRadius() ? p->getCutoffRadius() : maxCutoff;
            }
        }
        for (auto &&e : ctx->getAllOrder2Reactions()) {
            maxCutoff = maxCutoff < e->getEductDistance() ? e->getEductDistance() : maxCutoff;
        }
        NeighborList::maxCutoff = maxCutoff;
        if (maxCutoff > 0) {
            const auto desiredBoxWidth = .5 * maxCutoff;

            for (unsigned short i = 0; i < 3; ++i) {
                nBoxes[i] = static_cast<box_index>(floor(simBoxSize[i] / desiredBoxWidth));
                if (nBoxes[i] == 0) nBoxes[i] = 1;
                boxSize[i] = simBoxSize[i] / nBoxes[i];
            }
            for (box_index i = 0; i < nBoxes[0]; ++i) {
                for (box_index j = 0; j < nBoxes[1]; ++j) {
                    for (box_index k = 0; k < nBoxes[2]; ++k) {
                        boxes.push_back({i, j, k, nBoxes});
                    }
                }
            }
            for (box_index i = 0; i < nBoxes[0]; ++i) {
                for (box_index j = 0; j < nBoxes[1]; ++j) {
                    for (box_index k = 0; k < nBoxes[2]; ++k) {
                        setupNeighboringBoxes(i, j, k);
                    }
                }
            }
        }
    }
}

void NeighborList::setupNeighboringBoxes(const signed_box_index i, const signed_box_index j, const signed_box_index k) {
    auto me = getBox(i, j, k);
    for (signed_box_index _i = -2; _i < 3; ++_i) {
        for (signed_box_index _j = -2; _j < 3; ++_j) {
            for (signed_box_index _k = -2; _k < 3; ++_k) {
                // don't add me as neighbor to myself
                if (!(_i == 0 && _j == 0 && _k == 0)) {
                    me->addNeighbor(getBox(i + _i, j + _j, k + _k));
                }
            }
        }
    }
}

void NeighborList::clear() {
    boxes.clear();
}

void NeighborList::fillBoxes(const data_t &data) {
    if (maxCutoff > 0) {

        for (auto &&box : boxes) {
            box.particleIndices.clear();
        }
        pairs->clear();

        auto it_pos = data.cbegin_positions();
        particle_index idx = 0;
        while (it_pos != data.cend_positions()) {
            const auto pos = *it_pos;
            auto box = getBox(pos);
            if (box) {
                box->particleIndices.push_back(idx);
            }
            pairs->emplace(idx, std::vector<neighbor_t>());
            ++idx;
            ++it_pos;
        }
        {
            const auto size = boxes.size();
            const std::size_t grainSize = size / config->nThreads;
            const auto cutoffSquared = maxCutoff * maxCutoff;
            auto worker = [this, &data, cutoffSquared](const box_index begin, const box_index end) {
                const auto d2 = ctx->getDistSquaredFun();
                for (auto _b = begin; _b < end; ++_b) {
                    const auto &box = boxes[_b];
                    for (particle_index i = 0; i < box.particleIndices.size(); ++i) {
                        const auto pI = box.particleIndices[i];
                        if (pairs->find(pI) != pairs->end()) {
                            const auto pos = *(data.begin_positions() + pI);
                            for (particle_index j = 0; j < box.particleIndices.size(); ++j) {
                                if (i != j) {
                                    const auto pJ = box.particleIndices[j];
                                    const auto distSquared = d2(pos, *(data.begin_positions() + pJ));
                                    if (distSquared < cutoffSquared) {
                                        (*pairs)[pI].push_back({pJ, distSquared});
                                    }
                                }
                            }
                            for (auto &&neighboringBox : box.neighbors) {
                                for (const auto &pJ : neighboringBox->particleIndices) {
                                    const auto distSquared = d2(pos, *(data.begin_positions() + pJ));
                                    if (distSquared < cutoffSquared) {
                                        (*pairs)[pI].push_back({pJ, distSquared});
                                    }
                                }
                            }
                        } else {
                            log::console()->error("CPUNeighborList: The particle index was not to be found in the map, "
                                                          "risking a concurrent modification! This should not happen.");
                        }

                    }
                }
            };

            std::vector<util::ScopedThread> threads;
            threads.reserve(config->nThreads);

            for (auto i = 0; i < config->nThreads - 1; ++i) {
                threads.push_back(util::ScopedThread(std::thread(worker, i * grainSize, (i + 1) * grainSize)));
            }
            threads.push_back(
                    util::ScopedThread(std::thread(worker, (config->nThreads - 1) * grainSize, boxes.size())));
        }
    }
}

void NeighborList::create(const data_t &data) {
    simBoxSize = ctx->getBoxSize();
    setupBoxes();
    fillBoxes(data);
}

NeighborList::Box *NeighborList::getBox(signed_box_index i, signed_box_index j, signed_box_index k) {
    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nBoxes[0]);
    else if (i < 0 || i >= nBoxes[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nBoxes[1]);
    else if (j < 0 || j >= nBoxes[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nBoxes[2]);
    else if (k < 0 || k >= nBoxes[2]) return nullptr;
    return &boxes[k + j * nBoxes[2] + i * nBoxes[2] * nBoxes[1]];
}

void NeighborList::remove(const particle_index idx) {
    auto neighbors = (*pairs)[idx];
    for (auto &&neighbor : neighbors) {
        auto neighbors2 = (*pairs)[neighbor.idx];
        std::remove_if(neighbors2.begin(), neighbors2.end(), [idx](const neighbor_t &n) {
            return n.idx == idx;
        });
    }
    pairs->erase(idx);
}

void NeighborList::insert(const data_t &data, const particle_index idx) {
    const auto d2 = ctx->getDistSquaredFun();
    const auto pos = *(data.begin_positions() + idx);
    const auto cutoffSquared = maxCutoff * maxCutoff;
    auto box = getBox(pos);
    if (box) {
        box->particleIndices.push_back(idx);
        pairs->emplace(idx, std::vector<neighbor_t>());

        for (particle_index j = 0; j < box->particleIndices.size(); ++j) {
            if (idx != j) {
                const auto pJ = box->particleIndices[j];
                const auto distSquared = d2(pos, *(data.begin_positions() + pJ));
                if (distSquared < cutoffSquared) {
                    (*pairs)[idx].push_back({pJ, distSquared});
                    (*pairs)[pJ].push_back({idx, distSquared});
                }
            }
        }
        for (auto &&neighboringBox : box->neighbors) {
            for (const auto &pJ : neighboringBox->particleIndices) {
                const auto distSquared = d2(pos, *(data.begin_positions() + pJ));
                if (distSquared < cutoffSquared) {
                    (*pairs)[idx].push_back({pJ, distSquared});
                    (*pairs)[pJ].push_back({idx, distSquared});
                }
            }
        }
    } else {
        log::console()->error("could not assign particle (index={}) to any box!", idx);
    }
}

NeighborList::Box *NeighborList::getBox(const readdy::model::Particle::pos_type &pos) {
    const box_index i = static_cast<const box_index>(floor((pos[0] + .5 * simBoxSize[0]) / boxSize[0]));
    const box_index j = static_cast<const box_index>(floor((pos[1] + .5 * simBoxSize[1]) / boxSize[1]));
    const box_index k = static_cast<const box_index>(floor((pos[2] + .5 * simBoxSize[2]) / boxSize[2]));
    return getBox(i, j, k);
}

NeighborList::~NeighborList() = default;

Neighbor::Neighbor(const index_t idx, const double d2) : idx(idx), d2(d2) {}

}
}
}
}
