/**
 * << detailed description >>
 *
 * @file SingleCPUKernelStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 * @todo
 */

#include <readdy/model/Particle.h>
#include <vector>
#include <algorithm>
#include <readdy/kernel/singlecpu/SingleCPUKernelStateModel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            struct SingleCPUKernelStateModel::Impl {
                readdy::model::time_step_type t = 0;
                double currentEnergy = 0;
                bool firstUpdate = true;
                std::unique_ptr<model::SingleCPUParticleData> particleData;
                std::unique_ptr<model::SingleCPUNeighborList> neighborList;
                readdy::model::KernelContext const *context;
            };

            void SingleCPUKernelStateModel::updateModel(readdy::model::time_step_type t, bool forces) {
                const auto timeStepChanged = t != pimpl->t;
                const auto& difference = pimpl->context->getShortestDifferenceFun();
                pimpl->t = t;
                if (timeStepChanged || pimpl->firstUpdate) {
                    pimpl->neighborList->create(*pimpl->particleData);
                    pimpl->currentEnergy = 0;
                    pimpl->firstUpdate = false;
                    fireTimeStepChanged();
                }

                if (forces) {
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
                        for (auto &&it = pimpl->neighborList->begin(); it != pimpl->neighborList->end(); ++it) {
                            auto i = it->idx1;
                            auto j = it->idx2;
                            auto type_i = *(pimpl->particleData->begin_types() + i);
                            auto type_j = *(pimpl->particleData->begin_types() + j);
                            const auto &pos_i = *(pimpl->particleData->begin_positions() + i);
                            const auto &pos_j = *(pimpl->particleData->begin_positions() + j);
                            const auto &potentials = pimpl->context->getOrder2Potentials(type_i, type_j);
                            for (const auto &potential : potentials) {
                                potential->calculateForceAndEnergy(*(pimpl->particleData->begin_forces() + i), pimpl->currentEnergy, difference(pos_i, pos_j));
                            }
                        }
                    }
                }
            }


            SingleCPUKernelStateModel::SingleCPUKernelStateModel(readdy::model::KernelContext const *context) : pimpl(std::make_unique<SingleCPUKernelStateModel::Impl>()) {
                pimpl->particleData = std::make_unique<model::SingleCPUParticleData>(10);
                pimpl->neighborList = std::make_unique<model::NaiveSingleCPUNeighborList>();
                pimpl->context = context;
            }

            void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticle(const readdy::model::Particle &p) {
                pimpl->particleData->addParticles({p});
            }

            void readdy::kernel::singlecpu::SingleCPUKernelStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
                pimpl->particleData->addParticles(p);
            }

            const std::vector<readdy::model::Vec3> readdy::kernel::singlecpu::SingleCPUKernelStateModel::getParticlePositions() const {
                const auto size = pimpl->particleData->size();
                std::vector<readdy::model::Vec3> target{};
                target.reserve(size);
                std::copy(pimpl->particleData->begin_positions(), pimpl->particleData->begin_positions() + size, std::back_inserter(target));
                return target;
            }

            const readdy::model::time_step_type readdy::kernel::singlecpu::SingleCPUKernelStateModel::getCurrentTimeStep() const {
                return pimpl->t;
            }

            model::SingleCPUParticleData *readdy::kernel::singlecpu::SingleCPUKernelStateModel::getParticleData() const {
                return pimpl->particleData.get();
            }

            void readdy::kernel::singlecpu::SingleCPUKernelStateModel::removeParticle(const readdy::model::Particle &p) {
                pimpl->particleData->removeParticle(p);
            }

            double readdy::kernel::singlecpu::SingleCPUKernelStateModel::getEnergy() const {
                return pimpl->currentEnergy;
            }

            void readdy::kernel::singlecpu::SingleCPUKernelStateModel::increaseEnergy(double increase) {
                pimpl->currentEnergy += increase;
            }

            const model::SingleCPUNeighborList *const SingleCPUKernelStateModel::getNeighborList() const {
                return pimpl->neighborList.get();
            }

            const std::vector<readdy::model::Particle> SingleCPUKernelStateModel::getParticles() const {
                std::vector<readdy::model::Particle> result;
                for(auto i = 0; i < pimpl->particleData->size(); ++i) {
                    result.push_back((*pimpl->particleData)[i]);
                }
                return result;
            }


            SingleCPUKernelStateModel &SingleCPUKernelStateModel::operator=(SingleCPUKernelStateModel &&rhs) = default;

            SingleCPUKernelStateModel::SingleCPUKernelStateModel(SingleCPUKernelStateModel &&rhs) = default;

            SingleCPUKernelStateModel::~SingleCPUKernelStateModel() = default;
        }
    }
}




