/**
 * << detailed description >>
 *
 * @file CPUStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <readdy/kernel/cpu/CPUStateModel.h>
#include <future>

namespace readdy {
namespace kernel {
namespace cpu {

struct CPUStateModel::Impl {
    readdy::model::KernelContext *context;
    std::unique_ptr<readdy::kernel::cpu::model::NeighborList> neighborList;
    double currentEnergy = 0;

    template<bool fixpos=true>
    const model::ParticleData& cdata() const {
        if(fixpos) particleData->setFixPosFun(context->getFixPositionFun());
        return *particleData;
    }

    template<bool fixpos=true>
    model::ParticleData& data() {
        if(fixpos) particleData->setFixPosFun(context->getFixPositionFun());
        return *particleData;
    }

    Impl() {
        particleData = std::make_unique<CPUStateModel::data_t>();
    }

private:
    std::unique_ptr<readdy::kernel::cpu::model::ParticleData> particleData;
};

void CPUStateModel::calculateForces() {
    pimpl->currentEnergy = 0;
    const auto& particleData = pimpl->cdata<true>();
    // update forces and energy order 1 potentials
    {
        std::vector<double> energyUpdate;
        energyUpdate.reserve(config->nThreads);
        std::vector<util::scoped_thread> threads;
        threads.reserve(config->nThreads);
        using iter_t = data_t::entries_t::iterator;
        using pot1map = std::unordered_map<unsigned int, std::vector<readdy::model::potentials::PotentialOrder1 *>>;

        auto worker = [](iter_t begin, const iter_t end, double &energy,pot1map pot1Map) -> void {
            const readdy::model::Vec3 zero;
            for(auto it = begin; it != end; ++it) {
                if(!it->is_deactivated()) {
                    it->force = zero;
                    for (auto po1 : pot1Map[(*begin).type]) {
                        po1->calculateForceAndEnergy(it->force, energy, it->pos);
                    }
                }
            }
        };

        const auto size = particleData.size();
        const std::size_t grainSize = size / config->nThreads;

        auto it = particleData.begin();
        for (auto i = 0; i < config->nThreads - 1; ++i) {
            energyUpdate.push_back(0);
            threads.push_back(util::scoped_thread(
                    std::thread(worker, it, it + grainSize, std::ref(energyUpdate.back()),
                                pimpl->context->getAllOrder1Potentials())));
            it += grainSize;
        }
        energyUpdate.push_back(0);
        threads.push_back(
                util::scoped_thread(std::thread(worker, it, particleData.end(),
                                               std::ref(energyUpdate.back()), pimpl->context->getAllOrder1Potentials()))
        );
        std::for_each(energyUpdate.begin(), energyUpdate.end(), [this](double e) { pimpl->currentEnergy += e; });
    }

    // update forces and energy order 2 potentials
    {

        std::vector<std::future<double>> energyFutures;
        energyFutures.reserve(config->nThreads);
        {
            using iter_t = decltype(pimpl->neighborList->begin()->begin());
            using pot2map = std::unordered_map<readdy::util::ParticleTypePair, std::vector<readdy::model::potentials::PotentialOrder2 *>, readdy::util::ParticleTypePairHasher>;
            std::vector<util::scoped_thread> threads;
            threads.reserve(config->nThreads);
            const auto &d = pimpl->context->getShortestDifferenceFun();
            auto worker = [d](iter_t begin, iter_t end, std::promise<double> energyPromise, pot2map pot2Map) -> void {
                double energy = 0;

                for(auto it = begin; it != end; ++it) {
                    auto entry_i = it->first;

                    for (const auto &neighbor : it->second) {
                        readdy::model::Vec3 forceVec{0, 0, 0};
                        auto potit = pot2Map.find({entry_i->type, neighbor.idx->type});
                        if(potit != pot2Map.end()) {
                            for (const auto &potential : potit->second) {
                                if (neighbor.d2 < potential->getCutoffRadiusSquared()) {
                                    readdy::model::Vec3 updateVec{0, 0, 0};
                                    potential->calculateForceAndEnergy(updateVec, energy,
                                                                       d(entry_i->pos, neighbor.idx->pos));
                                    forceVec += updateVec;
                                }
                            }
                        }
                        entry_i->force += forceVec;
                    }
                }
                energyPromise.set_value(energy);
            };

            for (auto& map : *pimpl->neighborList) {
                std::promise<double> energyPromise;
                energyFutures.push_back(energyPromise.get_future());
                threads.push_back(util::scoped_thread(std::thread(worker, map.begin(), map.end(), std::move(energyPromise), pimpl->context->getAllOrder2Potentials())));
            }
        }
        for (auto &f : energyFutures) {
            pimpl->currentEnergy += f.get();
        }

    }
    /**
     * We need to take 0.5*energy as we iterate over every particle-neighbor-pair twice.
     */
    pimpl->currentEnergy /= 2.0;
}

const std::vector<readdy::model::Vec3> CPUStateModel::getParticlePositions() const {
    const auto& data = pimpl->cdata();
    std::vector<readdy::model::Vec3> target{};
    target.reserve(data.size());
    for(const auto& entry : data) {
        target.push_back(entry.position());
    }
    return target;
}

const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
    const auto& data = pimpl->cdata();
    std::vector<readdy::model::Particle> result;
    result.reserve(data.size());
    for(const auto& entry : data) {
        if(!entry.is_deactivated()) {
            result.push_back(data.toParticle(entry));
        }
    }
    return result;
}

void CPUStateModel::updateNeighborList() {
    pimpl->neighborList->create(pimpl->data<false>());
}

void CPUStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->data<false>().addParticle(p);
}

void CPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->data<false>().addParticles(p);
}

void CPUStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->data<false>().removeParticle(p);
}

double CPUStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

CPUStateModel::CPUStateModel(readdy::model::KernelContext *const context, util::Config const *const config)
        : pimpl(std::make_unique<Impl>()), config(config) {
    pimpl->context = context;
    pimpl->neighborList = std::make_unique<model::NeighborList>(context, config);
}

CPUStateModel::data_t *const CPUStateModel::getParticleData() const {
    return &pimpl->data<false>();
}

model::NeighborList *const CPUStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

void CPUStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void CPUStateModel::removeAllParticles() {
    pimpl->data<false>().clear();
}

CPUStateModel::~CPUStateModel() = default;


}
}
}