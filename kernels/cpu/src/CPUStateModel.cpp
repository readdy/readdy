/**
 * << detailed description >>
 *
 * @file CPUStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <readdy/kernel/singlecpu/model/SingleCPUParticleData.h>
#include <readdy/kernel/cpu/CPUStateModel.h>

namespace readdy {
    namespace kernel {
        namespace cpu {

            struct CPUStateModel::Impl {
                readdy::model::KernelContext *context;
                std::unique_ptr<readdy::kernel::singlecpu::model::SingleCPUParticleData> particleData;
                std::unique_ptr<readdy::kernel::cpu::model::NeighborList> neighborList;
                double currentEnergy = 0;
            };

            void CPUStateModel::calculateForces() {
                pimpl->currentEnergy = 0;
                // update forces and energy order 1 potentials
                std::vector<double> energyUpdate;
                energyUpdate.reserve(util::getNThreads());
                using data_t = readdy::kernel::singlecpu::model::SingleCPUParticleData;
                {
                    std::vector<util::ScopedThread> threads;
                    threads.reserve(util::getNThreads());
                    using _it_t = std::vector<unsigned int>::const_iterator;
                    using pot1map = std::unordered_map<unsigned int, std::vector<readdy::model::potentials::PotentialOrder1 *>>;

                    auto worker = [](_it_t it_types, _it_t it_types_end, double &energy,
                                     data_t *data, pot1map pot1Map) -> void {
                        auto idx = it_types - data->cbegin_types();
                        const readdy::model::Vec3 zero{0, 0, 0};
                        auto &&it_forces = data->begin_forces() + idx;
                        auto &&it_pos = data->cbegin_positions() + idx;
                        while (it_types != it_types_end) {
                            *it_forces = zero;
                            for (auto po1 : pot1Map[*it_types]) {
                                po1->calculateForceAndEnergy(*it_forces, energy, *it_pos);
                            }
                            ++it_forces;
                            ++it_types;
                            ++it_pos;
                        }
                    };

                    const auto size = pimpl->particleData->size();
                    const std::size_t grainSize = size / util::getNThreads();

                    auto it = pimpl->particleData->cbegin_types();
                    for (auto i = 0; i < util::getNThreads() - 1; ++i) {
                        energyUpdate.push_back(0);
                        threads.push_back(util::ScopedThread(
                                std::thread(worker, it, it + grainSize, std::ref(energyUpdate.back()),
                                            pimpl->particleData.get(), pimpl->context->getAllOrder1Potentials())));
                        it += grainSize;
                    }
                    energyUpdate.push_back(0);
                    threads.push_back(
                            util::ScopedThread(std::thread(worker, it, pimpl->particleData->cend_types(),
                                                           std::ref(energyUpdate.back()), pimpl->particleData.get(),
                                                           pimpl->context->getAllOrder1Potentials()))
                    );
                }


                // update forces and energy order 2 potentials
                {
                    using nl_it_t = model::NeighborList::container_t::iterator;
                    using pot2map = std::unordered_map<readdy::util::ParticleTypePair, std::vector<readdy::model::potentials::PotentialOrder2*>, readdy::util::ParticleTypePairHasher>;
                    using dist_t = std::function<readdy::model::Vec3(const readdy::model::Vec3 &, const readdy::model::Vec3 &)>;
                    std::vector<util::ScopedThread> threads;
                    threads.reserve(util::getNThreads());
                    auto worker = [](nl_it_t begin, const unsigned long n, double& energy, data_t *data, pot2map pot2Map, dist_t dist) -> void {

                        auto it = begin;
                        for (unsigned long _i = 0; _i < n; ++_i) {
                            auto i = it->first;

                            for (auto &&j : it->second) {
                                readdy::model::Vec3 forceVec{0, 0, 0};
                                const auto type_i = *(data->cbegin_types() + i);
                                const auto type_j = *(data->cbegin_types() + j);
                                const auto &pos_i = *(data->cbegin_positions() + i);
                                const auto &pos_j = *(data->cbegin_positions() + j);
                                const auto &potentials = pot2Map[{type_i, type_j}];
                                for (const auto &potential : potentials) {
                                    readdy::model::Vec3 updateVec{0, 0, 0};
                                    potential->calculateForceAndEnergy(updateVec, energy,
                                                                       dist(pos_i, pos_j));
                                    forceVec += updateVec;
                                }
                                *(data->begin_forces() + i) += forceVec;
                            }

                            it = std::next(it);
                        }
                    };

                    const auto size = pimpl->neighborList->pairs->size();
                    const std::size_t grainSize = size / util::getNThreads();

                    auto it = pimpl->neighborList->pairs->begin();
                    for (auto i = 0; i < util::getNThreads() - 1; ++i) {
                        threads.push_back(util::ScopedThread(std::thread(worker, it, grainSize, std::ref(energyUpdate[i]), pimpl->particleData.get(), pimpl->context->getAllOrder2Potentials(), pimpl->context->getShortestDifferenceFun())));
                        std::advance(it, grainSize);
                    }
                    threads.push_back(util::ScopedThread(
                            std::thread(worker, it, std::distance(it, pimpl->neighborList->pairs->end()), std::ref(energyUpdate.back()), pimpl->particleData.get(), pimpl->context->getAllOrder2Potentials(), pimpl->context->getShortestDifferenceFun())));

                }
                std::for_each(energyUpdate.begin(), energyUpdate.end(),
                              [this](const double &e) { pimpl->currentEnergy += e; });
                /**
                 * We need to take 0.5*energy as we iterate over every particle-neighbor-pair twice.
                 */
                pimpl->currentEnergy /= 2.0;
            }

            const std::vector<readdy::model::Vec3> CPUStateModel::getParticlePositions() const {
                const auto begin = pimpl->particleData->begin_positions();
                const auto end = pimpl->particleData->end_positions();
                std::vector<readdy::model::Vec3> target{};
                target.reserve(end - begin);
                std::copy(begin, end, std::back_inserter(target));
                return target;
            }

            const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
                std::vector<readdy::model::Particle> result;
                for (auto i = 0; i < pimpl->particleData->size(); ++i) {
                    result.push_back((*pimpl->particleData)[i]);
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

            CPUStateModel::CPUStateModel(readdy::model::KernelContext *const context) : pimpl(
                    std::make_unique<Impl>()) {
                pimpl->context = context;
                pimpl->particleData = std::make_unique<readdy::kernel::singlecpu::model::SingleCPUParticleData>(0);
                pimpl->neighborList = std::make_unique<model::NeighborList>(context);
            }

            readdy::kernel::singlecpu::model::SingleCPUParticleData *const CPUStateModel::getParticleData() const {
                return pimpl->particleData.get();
            }

            const model::NeighborList *const CPUStateModel::getNeighborList() const {
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