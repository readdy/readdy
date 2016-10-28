/**
 * << detailed description >>
 *
 * @file CPUStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <readdy/kernel/cpu/CPUStateModel.h>
#include <readdy/kernel/cpu/model/ParticleData.h>

namespace readdy {
namespace kernel {
namespace cpu {

struct CPUStateModel::Impl {
    readdy::model::KernelContext *context;
    std::unique_ptr<readdy::kernel::cpu::model::ParticleData> particleData;
    std::unique_ptr<readdy::kernel::cpu::model::NeighborList> neighborList;
    double currentEnergy = 0;
};

void CPUStateModel::calculateForces() {
    pimpl->currentEnergy = 0;
    // update forces and energy order 1 potentials
    std::vector<double> energyUpdate;
    energyUpdate.reserve(config->nThreads);
    {
        std::vector<util::scoped_thread> threads;
        threads.reserve(config->nThreads);
        using iter_t = data_t::entries_t::iterator;
        using pot1map = std::unordered_map<unsigned int, std::vector<readdy::model::potentials::PotentialOrder1 *>>;

        auto worker = [](iter_t begin, const iter_t end, double &energy,pot1map pot1Map) -> void {
            const readdy::model::Vec3 zero;
            while (begin != end) {
                if(!begin->is_deactivated()) {
                    (*begin).force = zero;
                    for (auto po1 : pot1Map[(*begin).type]) {
                        po1->calculateForceAndEnergy((*begin).force, energy, (*begin).pos);
                    }
                }
                ++begin;
            }
        };

        const auto size = pimpl->particleData->size();
        const std::size_t grainSize = size / config->nThreads;

        auto it = pimpl->particleData->entries.begin();
        for (auto i = 0; i < config->nThreads - 1; ++i) {
            energyUpdate.push_back(0);
            threads.push_back(util::scoped_thread(
                    std::thread(worker, it, it + grainSize, std::ref(energyUpdate.back()),
                                pimpl->context->getAllOrder1Potentials())));
            it += grainSize;
        }
        energyUpdate.push_back(0);
        threads.push_back(
                util::scoped_thread(std::thread(worker, it, pimpl->particleData->entries.end(),
                                               std::ref(energyUpdate.back()), pimpl->context->getAllOrder1Potentials()))
        );
    }


    // update forces and energy order 2 potentials
    {
        using nl_it_t = model::NeighborList::container_t::iterator;
        using pot2map = std::unordered_map<readdy::util::ParticleTypePair, std::vector<readdy::model::potentials::PotentialOrder2 *>, readdy::util::ParticleTypePairHasher>;
        using dist_t = std::function<readdy::model::Vec3(const readdy::model::Vec3 &, const readdy::model::Vec3 &)>;
        std::vector<util::scoped_thread> threads;
        threads.reserve(config->nThreads);
        auto worker = [](nl_it_t begin, const unsigned long n, double &energy, data_t *data, pot2map pot2Map,
                         dist_t dist) -> void {
            auto it = begin;
            for (unsigned long _i = 0; _i < n; ++_i) {
                auto entry_i = it->first;

                for (const auto &neighbor : it->second) {
                    readdy::model::Vec3 forceVec{0, 0, 0};
                    const auto &potentials = pot2Map[{entry_i->type, neighbor.idx->type}];
                    for (const auto &potential : potentials) {
                        if (neighbor.d2 < potential->getCutoffRadiusSquared()) {
                            readdy::model::Vec3 updateVec{0, 0, 0};
                            potential->calculateForceAndEnergy(updateVec, energy, dist(entry_i->pos, neighbor.idx->pos));
                            forceVec += updateVec;
                        }
                    }
                    entry_i->force += forceVec;
                }

                it = std::next(it);
            }
        };

        const auto size = pimpl->neighborList->pairs.size();
        const std::size_t grainSize = size / config->nThreads;

        auto it = pimpl->neighborList->pairs.begin();
        for (auto i = 0; i < config->nThreads - 1; ++i) {
            threads.push_back(util::scoped_thread(
                    std::thread(worker, it, grainSize, std::ref(energyUpdate[i]), pimpl->particleData.get(),
                                pimpl->context->getAllOrder2Potentials(), pimpl->context->getShortestDifferenceFun())));
            std::advance(it, grainSize);
        }
        threads.push_back(util::scoped_thread(
                std::thread(worker, it, std::distance(it, pimpl->neighborList->pairs.end()),
                            std::ref(energyUpdate.back()), pimpl->particleData.get(),
                            pimpl->context->getAllOrder2Potentials(), pimpl->context->getShortestDifferenceFun())));

    }
    std::for_each(energyUpdate.begin(), energyUpdate.end(),
                  [this](const double &e) { pimpl->currentEnergy += e; });
    /**
     * We need to take 0.5*energy as we iterate over every particle-neighbor-pair twice.
     */
    pimpl->currentEnergy /= 2.0;
}

const std::vector<readdy::model::Vec3> CPUStateModel::getParticlePositions() const {
    const auto& entries = pimpl->particleData->entries;
    std::vector<readdy::model::Vec3> target{};
    target.reserve(pimpl->particleData->size());
    std::for_each(entries.begin(), entries.end(), [&target](const data_t::Entry& entry) {
        if(!entry.is_deactivated()) {
            target.push_back(entry.pos);
        }
    });
    return target;
}

const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
    std::vector<readdy::model::Particle> result;
    result.reserve(pimpl->particleData->size());
    for(const auto& entry : pimpl->particleData->entries) {
        if(!entry.is_deactivated()) {
            result.push_back(pimpl->particleData->toParticle(entry));
        }
    }
    return result;
}

void CPUStateModel::updateNeighborList() {
    pimpl->neighborList->create(*pimpl->particleData);
}

void CPUStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->particleData->addParticle(p);
}

void CPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->particleData->addParticles(p);
}

void CPUStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->particleData->removeParticle(p);
}

double CPUStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

CPUStateModel::CPUStateModel(readdy::model::KernelContext *const context, util::Config const *const config)
        : pimpl(std::make_unique<Impl>()), config(config) {
    pimpl->context = context;
    pimpl->particleData = std::make_unique<data_t>();
    pimpl->neighborList = std::make_unique<model::NeighborList>(context, config);
}

CPUStateModel::data_t *const CPUStateModel::getParticleData() const {
    return pimpl->particleData.get();
}

model::NeighborList *const CPUStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

void CPUStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void CPUStateModel::removeAllParticles() {
    pimpl->particleData->clear();
}

CPUStateModel::~CPUStateModel() = default;


}
}
}