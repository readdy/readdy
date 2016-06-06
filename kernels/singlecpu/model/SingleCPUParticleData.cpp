/**
 * << detailed description >>
 *
 * @file SingleCPUParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#include <readdy/common/make_unique.h>
#include <numeric>
#include <boost/log/trivial.hpp>
#include "SingleCPUParticleData.h"

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace model {
                SingleCPUParticleData::SingleCPUParticleData() : SingleCPUParticleData(0) {
                }

                SingleCPUParticleData::SingleCPUParticleData(unsigned int capacity) {
                    ids = std::make_unique<std::vector<boost::uuids::uuid>>(capacity);
                    positions = std::make_unique<std::vector<readdy::model::Vec3>>(capacity);
                    forces = std::make_unique<std::vector<readdy::model::Vec3>>(capacity);
                    type = std::make_unique<std::vector<unsigned int>>(capacity);
                    deactivatedParticles = std::make_unique<std::vector<size_t>>(capacity);
                    std::iota(deactivatedParticles->begin(), deactivatedParticles->end(), 0);
                }

                void SingleCPUParticleData::swap(SingleCPUParticleData &rhs) {
                    std::swap(ids, rhs.ids);
                    std::swap(positions, rhs.positions);
                    std::swap(forces, rhs.forces);
                    std::swap(type, rhs.type);
                    std::swap(deactivatedParticles, rhs.deactivatedParticles);
                }

                size_t SingleCPUParticleData::size() {
                    return ids->size() - deactivatedParticles->size();
                }

                size_t SingleCPUParticleData::max_size() {
                    return ids->max_size();
                }

                bool SingleCPUParticleData::empty() {
                    return size() == 0;
                }

                void SingleCPUParticleData::addParticles(const std::vector<readdy::model::Particle> &particles) {
                    auto added = particles.cbegin();
                    while (added != particles.cend()) {
                        if (deactivatedParticles->size() > 0) {
                            const auto idx = (*deactivatedParticles)[0];

                            (*ids)[idx] = added->getId();
                            (*positions)[idx] = added->getPos();
                            (*forces)[idx] = {0, 0, 0};
                            (*type)[idx] = added->getType();

                            deactivatedParticles->erase(deactivatedParticles->begin());
                        } else {
                            ids->push_back(added->getId());
                            positions->push_back(added->getPos());
                            forces->push_back({0, 0, 0});
                            type->push_back(added->getType());
                        }
                        added++;
                    }
                }

                void SingleCPUParticleData::removeParticle(const size_t index) {
                    const auto pos = (begin_ids() + index).getInternalPosition();
                    deactivatedParticles->push_back(pos);
                }

                void SingleCPUParticleData::removeParticle(const readdy::model::Particle &particle) {
                    auto &&beginIt = begin_ids();
                    auto &&endIt = end_ids();
                    auto &&it = std::find(beginIt, endIt, particle.getId());
                    if (it != endIt) {
                        deactivatedParticles->push_back(it.getInternalPosition());
                    } else {
                        BOOST_LOG_TRIVIAL(warning) << "Could not find and thus remove particle";
                    }
                }

                SingleCPUParticleData::skipping_iterator<boost::uuids::uuid> SingleCPUParticleData::begin_ids() {
                    return SingleCPUParticleData::skipping_iterator<boost::uuids::uuid>(this, 0, ids->begin());
                }

                SingleCPUParticleData::skipping_iterator<boost::uuids::uuid> SingleCPUParticleData::end_ids() {
                    return SingleCPUParticleData::skipping_iterator<boost::uuids::uuid>(this, ids->size(), ids->end());
                }

                SingleCPUParticleData::skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::begin_positions() {
                    return SingleCPUParticleData::skipping_iterator<readdy::model::Vec3>(this, 0, positions->begin());
                }

                SingleCPUParticleData::skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::end_positions() {
                    return SingleCPUParticleData::skipping_iterator<readdy::model::Vec3>(this, positions->size(), positions->end());
                }

                SingleCPUParticleData::skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::begin_forces() {
                    return SingleCPUParticleData::skipping_iterator<readdy::model::Vec3>(this, 0, forces->begin());
                }

                SingleCPUParticleData::skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::end_forces() {
                    return SingleCPUParticleData::skipping_iterator<readdy::model::Vec3>(this, forces->size(), forces->end());
                }

                SingleCPUParticleData::skipping_iterator<unsigned int> SingleCPUParticleData::begin_types() {
                    return SingleCPUParticleData::skipping_iterator<unsigned int>(this, 0, type->begin());
                }

                SingleCPUParticleData::skipping_iterator<unsigned int> SingleCPUParticleData::end_types() {
                    return SingleCPUParticleData::skipping_iterator<unsigned int>(this, type->size(), type->end());
                }

                readdy::model::Particle SingleCPUParticleData::operator[](const size_t index) {
                    return readdy::model::Particle(*(begin_positions() + index), *(begin_types() + index), *(begin_ids() + index));
                }


                SingleCPUParticleData &SingleCPUParticleData::operator=(SingleCPUParticleData &&rhs) = default;

                SingleCPUParticleData::SingleCPUParticleData(SingleCPUParticleData &&rhs) = default;

                SingleCPUParticleData::~SingleCPUParticleData() {
                }

                std::vector<size_t> *SingleCPUParticleData::getDeactivatedParticles() const {
                    return deactivatedParticles.get();
                };


            }
        }
    }
}




