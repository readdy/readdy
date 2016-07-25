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
#include <readdy/kernel/cpu/model/NeighborList.h>

namespace readdy {
    namespace kernel {
        namespace cpu {

            struct CPUStateModel::Impl {
                readdy::model::KernelContext* context;
                std::unique_ptr<readdy::kernel::singlecpu::model::SingleCPUParticleData> particleData;
                std::unique_ptr<readdy::kernel::cpu::model::NeighborList> neighborList;
                double currentEnergy = 0;
            };

            void CPUStateModel::calculateForces() {
                // update forces and energy order 1 potentials
                {
                    const readdy::model::Vec3 zero{0, 0, 0};
                    auto&& it_forces = pimpl->particleData->begin_forces();
                    auto&& it_types = pimpl->particleData->begin_types();
                    auto&& it_pos = pimpl->particleData->begin_positions();
                    while(it_forces != pimpl->particleData->end_forces()) {
                        *it_forces = zero;
                        for(const auto &po1 : pimpl->context->getOrder1Potentials(*it_types)) {
                            po1->calculateForceAndEnergy(*it_forces, pimpl->currentEnergy, *it_pos);
                        }
                        ++it_forces;
                        ++it_types;
                        ++it_pos;
                    }
                }

                // update forces and energy order 2 potentials
                {
                    const auto& difference = pimpl->context->getShortestDifferenceFun();
                    readdy::model::Vec3 forceVec {0,0,0};
                    for (auto &&it = pimpl->neighborList->begin(); it != pimpl->neighborList->end(); ++it) {
                        auto i = it->idx1;
                        auto j = it->idx2;
                        auto type_i = *(pimpl->particleData->begin_types() + i);
                        auto type_j = *(pimpl->particleData->begin_types() + j);
                        const auto &pos_i = *(pimpl->particleData->begin_positions() + i);
                        const auto &pos_j = *(pimpl->particleData->begin_positions() + j);
                        const auto &potentials = pimpl->context->getOrder2Potentials(type_i, type_j);
                        for (const auto &potential : potentials) {
                            potential->calculateForceAndEnergy(forceVec, pimpl->currentEnergy, difference(pos_i, pos_j));
                            *(pimpl->particleData->begin_forces() + i) += forceVec;
                            *(pimpl->particleData->begin_forces() + j) += -1*forceVec;
                        }
                    }
                }
            }

            const std::vector<readdy::model::Vec3> CPUStateModel::getParticlePositions() const {
                const auto begin = pimpl->particleData->begin_positions();
                const auto end = pimpl->particleData->end_positions();
                std::vector<readdy::model::Vec3> target{};
                target.reserve(end-begin);
                std::copy(begin, end, std::back_inserter(target));
                return target;
            }

            const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
                std::vector<readdy::model::Particle> result;
                for(auto i = 0; i < pimpl->particleData->size(); ++i) {
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

            CPUStateModel::CPUStateModel(readdy::model::KernelContext *const context) : pimpl(std::make_unique<Impl>())
            {
                pimpl->context = context;
                pimpl->particleData = std::make_unique<readdy::kernel::singlecpu::model::SingleCPUParticleData>(0);
                pimpl->neighborList = std::make_unique<model::NeighborList>(context);
            }

            std::vector<model::ParticleIndexPair>::iterator CPUStateModel::begin_neighborList() {
                return pimpl->neighborList->begin();
            }

            std::vector<model::ParticleIndexPair>::iterator CPUStateModel::end_neighborList() {
                return pimpl->neighborList->end();
            }

            std::vector<model::ParticleIndexPair>::const_iterator CPUStateModel::begin_neighborList() const {
                return cbegin_neighborList();
            }

            std::vector<model::ParticleIndexPair>::const_iterator CPUStateModel::end_neighborList() const {
                return cend_neighborList();
            }

            std::vector<model::ParticleIndexPair>::const_iterator CPUStateModel::cbegin_neighborList() const {
                return pimpl->neighborList->cbegin();
            }

            std::vector<model::ParticleIndexPair>::const_iterator CPUStateModel::cend_neighborList() const {
                return pimpl->neighborList->cend();
            }

            readdy::kernel::singlecpu::model::SingleCPUParticleData *const CPUStateModel::getParticleData() const {
                return pimpl->particleData.get();
            }


            CPUStateModel::~CPUStateModel() = default;


        }
    }
}