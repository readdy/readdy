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
#include <readdy/kernel/singlecpu/model/SingleCPUParticleData.h>

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
                    deactivatedParticles = std::make_unique<boost::container::flat_set<size_t>>();
                    auto&& hint = deactivatedParticles->begin();
                    for(size_t i = 0; i < capacity; i++) {
                        hint = deactivatedParticles->insert(hint, i);
                    }
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

                void SingleCPUParticleData::addParticle(const readdy::model::Particle &particle) {
                    addParticles({particle});
                };

                void SingleCPUParticleData::addParticles(const std::vector<readdy::model::Particle> &particles) {
                    auto added = particles.cbegin();
                    while (added != particles.cend()) {
                        if (deactivatedParticles->size() > 0) {
                            const auto begin = deactivatedParticles->begin();
                            const auto idx = *begin;

                            (*ids)[idx] = added->getId();
                            (*positions)[idx] = added->getPos();
                            (*forces)[idx] = {0, 0, 0};
                            (*type)[idx] = added->getType();

                            deactivatedParticles->erase(begin);
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
                    deactivatedParticles->insert(pos);
                }

                void SingleCPUParticleData::removeParticle(const readdy::model::Particle &particle) {
                    auto &&beginIt = begin_ids();
                    auto &&endIt = end_ids();
                    auto &&it = std::find(beginIt, endIt, particle.getId());
                    if (it != endIt) {
                        deactivatedParticles->insert(it.getInternalPosition());
                    } else {
                        BOOST_LOG_TRIVIAL(warning) << "Could not find and thus remove particle";
                    }
                }

                SingleCPUParticleData::skipping_iterator<boost::uuids::uuid> SingleCPUParticleData::begin_ids() {
                    return SingleCPUParticleData::skipping_iterator<boost::uuids::uuid>(this, ids->begin(), ids->begin());
                }

                SingleCPUParticleData::skipping_iterator<boost::uuids::uuid> SingleCPUParticleData::end_ids() {
                    return SingleCPUParticleData::skipping_iterator<boost::uuids::uuid>(this, ids->begin(), ids->end());
                }

                SingleCPUParticleData::skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::begin_positions() {
                    return SingleCPUParticleData::skipping_iterator<readdy::model::Vec3>(this, positions->begin(), positions->begin());
                }

                SingleCPUParticleData::skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::end_positions() {
                    return SingleCPUParticleData::skipping_iterator<readdy::model::Vec3>(this, positions->begin(), positions->end());
                }

                SingleCPUParticleData::skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::begin_forces() {
                    return SingleCPUParticleData::skipping_iterator<readdy::model::Vec3>(this, forces->begin(), forces->begin());
                }

                SingleCPUParticleData::skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::end_forces() {
                    return SingleCPUParticleData::skipping_iterator<readdy::model::Vec3>(this, forces->begin(), forces->end());
                }

                SingleCPUParticleData::skipping_iterator<unsigned int> SingleCPUParticleData::begin_types() {
                    return SingleCPUParticleData::skipping_iterator<unsigned int>(this, type->begin(), type->begin());
                }

                SingleCPUParticleData::skipping_iterator<unsigned int> SingleCPUParticleData::end_types() {
                    return SingleCPUParticleData::skipping_iterator<unsigned int>(this, type->begin(), type->end());
                }

                readdy::model::Particle SingleCPUParticleData::operator[](const size_t index) {
                    return readdy::model::Particle(*(begin_positions() + index), *(begin_types() + index), *(begin_ids() + index));
                }


                SingleCPUParticleData &SingleCPUParticleData::operator=(SingleCPUParticleData &&rhs) = default;

                SingleCPUParticleData::SingleCPUParticleData(SingleCPUParticleData &&rhs) = default;

                SingleCPUParticleData::~SingleCPUParticleData() {
                }

                boost::container::flat_set<size_t> *SingleCPUParticleData::getDeactivatedParticles() const {
                    return deactivatedParticles.get();
                }

                SingleCPUParticleData::const_skipping_iterator <boost::uuids::uuid> SingleCPUParticleData::begin_ids() const {
                    return cbegin_ids();
                }

                SingleCPUParticleData::const_skipping_iterator<boost::uuids::uuid> SingleCPUParticleData::cbegin_ids() const {
                    return SingleCPUParticleData::const_skipping_iterator<boost::uuids::uuid>(this, ids->cbegin(), ids->cbegin());
                }

                SingleCPUParticleData::const_skipping_iterator<boost::uuids::uuid> SingleCPUParticleData::end_ids() const {
                    return cend_ids();
                }

                SingleCPUParticleData::const_skipping_iterator<boost::uuids::uuid> SingleCPUParticleData::cend_ids() const {
                    return SingleCPUParticleData::const_skipping_iterator<boost::uuids::uuid>(this, ids->cbegin(), ids->cend());
                }

                SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::begin_positions() const {
                    return cbegin_positions();
                }

                SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::cbegin_positions() const {
                    return SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3>(this, positions->cbegin(), positions->cbegin());
                }

                SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::end_positions() const {
                    return cend_positions();
                }

                SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::cend_positions() const {
                    return SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3>(this, positions->cbegin(), positions->cend());
                }

                SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::begin_forces() const {
                    return cbegin_forces();
                }

                SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::cbegin_forces() const {
                    return const_skipping_iterator<readdy::model::Vec3>(this, forces->cbegin(), forces->cbegin());
                }

                SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::end_forces() const {
                    return cend_forces();
                }

                SingleCPUParticleData::const_skipping_iterator<readdy::model::Vec3> SingleCPUParticleData::cend_forces() const {
                    return const_skipping_iterator<readdy::model::Vec3>(this, forces->cbegin(), forces->cend());
                }

                SingleCPUParticleData::const_skipping_iterator<unsigned int> SingleCPUParticleData::begin_types() const {
                    return cbegin_types();
                }

                SingleCPUParticleData::const_skipping_iterator<unsigned int> SingleCPUParticleData::cbegin_types() const {
                    return SingleCPUParticleData::const_skipping_iterator<unsigned int>(this, type->cbegin(), type->cbegin());
                }

                SingleCPUParticleData::const_skipping_iterator<unsigned int> SingleCPUParticleData::end_types() const {
                    return cend_types();
                }

                SingleCPUParticleData::const_skipping_iterator<unsigned int> SingleCPUParticleData::cend_types() const {
                    return SingleCPUParticleData::const_skipping_iterator<unsigned int>(this, type->cbegin(), type->cend());
                }


            }
        }
    }
}




