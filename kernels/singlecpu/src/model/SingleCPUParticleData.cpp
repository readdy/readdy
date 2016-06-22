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
                    markedForDeactivation = std::make_unique<std::set<size_t>>();
                    deactivated = std::make_unique<std::vector<bool>>(capacity);
                    std::fill(deactivated->begin(), deactivated->end(), true);
                    n_deactivated = capacity;
                    deactivated_index = 0;
                }

                void SingleCPUParticleData::swap(SingleCPUParticleData &rhs) {
                    std::swap(ids, rhs.ids);
                    std::swap(positions, rhs.positions);
                    std::swap(forces, rhs.forces);
                    std::swap(type, rhs.type);
                    std::swap(deactivated, rhs.deactivated);
                    std::swap(deactivated_index, rhs.deactivated_index);
                    std::swap(n_deactivated, rhs.n_deactivated);
                    std::swap(markedForDeactivation, rhs.markedForDeactivation);
                }

                size_t SingleCPUParticleData::size() const {
                    auto s = markedForDeactivation->size();
                    return s <= deactivated_index ? deactivated_index - s : 0;
                }

                size_t SingleCPUParticleData::max_size() const {
                    return ids->max_size();
                }

                bool SingleCPUParticleData::empty() const {
                    return size() == 0;
                }

                void SingleCPUParticleData::addParticle(const readdy::model::Particle &particle) {
                    addParticles({particle});
                };

                void SingleCPUParticleData::addParticles(const std::vector<readdy::model::Particle> &particles) {
                    auto added = particles.cbegin();
                    auto ids_it = ids->begin() + deactivated_index;
                    auto positions_it = positions->begin() + deactivated_index;
                    auto forces_it = forces->begin() + deactivated_index;
                    auto type_it = type->begin() + deactivated_index;
                    auto deactivated_it = deactivated->begin() + deactivated_index;
                    while (added != particles.cend()) {
                        if (n_deactivated > 0) {

                            *ids_it = added->getId();
                            *positions_it = added->getPos();
                            *forces_it = {0, 0, 0};
                            *type_it = added->getType();
                            *deactivated_it = false;

                            --n_deactivated;
                            ++deactivated_index;

                            ++ids_it; ++positions_it; ++forces_it; ++type_it; ++deactivated_it;
                        } else {
                            ids->push_back(added->getId());
                            positions->push_back(added->getPos());
                            forces->push_back({0, 0, 0});
                            type->push_back(added->getType());
                            deactivated->push_back(false);
                            ++deactivated_index;
                        }
                        ++added;
                    }
                }

                void SingleCPUParticleData::removeParticle(const size_t index) {
                    (*deactivated)[index] = true;

                    std::swap((*ids)[index], (*ids)[deactivated_index-1]);
                    std::swap((*positions)[index], (*positions)[deactivated_index-1]);
                    std::swap((*forces)[index], (*forces)[deactivated_index-1]);
                    std::swap((*type)[index], (*type)[deactivated_index-1]);
                    std::swap((*deactivated)[index], (*deactivated)[deactivated_index-1]);

                    ++n_deactivated;
                    if(deactivated_index == 0) throw std::runtime_error("hier sollte man aber nicht hinkommen!1");
                    --deactivated_index;
                }

                void SingleCPUParticleData::markForDeactivation(size_t index) {
                    (*deactivated)[index] = true;
                    markedForDeactivation->emplace(index);
                }

                void SingleCPUParticleData::deactivateMarked() {
                    if(markedForDeactivation->size() == 0) return;
                    if(deactivated->size() < deactivated_index-1) {
                        throw std::runtime_error("this should not happen");
                    }
                    if (deactivated_index > 0) {
                        auto deactivatedIt = deactivated->begin() + deactivated_index - 1;
                        for (auto &&idx : *markedForDeactivation) {
                            while (*deactivatedIt && deactivatedIt != deactivated->begin()) {
                                --deactivated_index;
                                ++n_deactivated;
                                --deactivatedIt;
                            }
                            if (idx < deactivated_index) {
                                removeParticle(idx);
                                --deactivatedIt;
                            } else {
                                break;
                            }
                        }
                    }
                    markedForDeactivation->clear();
                }

                void SingleCPUParticleData::removeParticle(const readdy::model::Particle &particle) {
                    auto &&beginIt = begin_ids();
                    auto &&endIt = end_ids();
                    auto &&it = std::find(beginIt, endIt, particle.getId());
                    if (it != endIt) {
                        removeParticle(it-beginIt);
                    } else {
                        BOOST_LOG_TRIVIAL(warning) << "Could not find and thus remove particle";
                    }
                }

                std::vector<boost::uuids::uuid>::iterator SingleCPUParticleData::begin_ids() {
                    return ids->begin();
                }

                std::vector<boost::uuids::uuid>::iterator SingleCPUParticleData::end_ids() {
                    return ids->begin() + deactivated_index;
                }

                std::vector<readdy::model::Vec3>::iterator SingleCPUParticleData::begin_positions() {
                    return positions->begin();
                }

                std::vector<readdy::model::Vec3>::iterator SingleCPUParticleData::end_positions() {
                    return positions->begin()+deactivated_index;
                }

                std::vector<readdy::model::Vec3>::iterator SingleCPUParticleData::begin_forces() {
                    return forces->begin();
                }

                std::vector<readdy::model::Vec3>::iterator SingleCPUParticleData::end_forces() {
                    return forces->begin()+deactivated_index;
                }

                std::vector<unsigned int>::iterator SingleCPUParticleData::begin_types() {
                    return type->begin();
                }

                std::vector<unsigned int>::iterator SingleCPUParticleData::end_types() {
                    return type->begin()+deactivated_index;
                }

                readdy::model::Particle SingleCPUParticleData::operator[](const size_t index) const{
                    return readdy::model::Particle(*(begin_positions() + index), *(begin_types() + index), *(begin_ids() + index));
                }


                SingleCPUParticleData &SingleCPUParticleData::operator=(SingleCPUParticleData &&rhs) = default;

                SingleCPUParticleData::SingleCPUParticleData(SingleCPUParticleData &&rhs) = default;

                SingleCPUParticleData::~SingleCPUParticleData() {
                }

                std::vector <boost::uuids::uuid>::const_iterator SingleCPUParticleData::begin_ids() const {
                    return cbegin_ids();
                }

                std::vector<boost::uuids::uuid>::const_iterator SingleCPUParticleData::cbegin_ids() const {
                    return ids->cbegin();
                }

                std::vector<boost::uuids::uuid>::const_iterator SingleCPUParticleData::end_ids() const {
                    return cend_ids();
                }

                std::vector<boost::uuids::uuid>::const_iterator SingleCPUParticleData::cend_ids() const {
                    return ids->begin()+deactivated_index;
                }

                std::vector<readdy::model::Vec3>::const_iterator SingleCPUParticleData::begin_positions() const {
                    return cbegin_positions();
                }

                std::vector<readdy::model::Vec3>::const_iterator SingleCPUParticleData::cbegin_positions() const {
                    return positions->cbegin();
                }

                std::vector<readdy::model::Vec3>::const_iterator SingleCPUParticleData::end_positions() const {
                    return cend_positions();
                }

                std::vector<readdy::model::Vec3>::const_iterator SingleCPUParticleData::cend_positions() const {
                   return positions->cbegin()+deactivated_index;
                }

                std::vector<readdy::model::Vec3>::const_iterator SingleCPUParticleData::begin_forces() const {
                    return cbegin_forces();
                }

                std::vector<readdy::model::Vec3>::const_iterator SingleCPUParticleData::cbegin_forces() const {
                    return forces->cbegin();
                }

                std::vector<readdy::model::Vec3>::const_iterator SingleCPUParticleData::end_forces() const {
                    return cend_forces();
                }

                std::vector<readdy::model::Vec3>::const_iterator SingleCPUParticleData::cend_forces() const {
                    return forces->cbegin()+deactivated_index;
                }

                std::vector<unsigned int>::const_iterator SingleCPUParticleData::begin_types() const {
                    return cbegin_types();
                }

                std::vector<unsigned int>::const_iterator SingleCPUParticleData::cbegin_types() const {
                    return type->cbegin();
                }

                std::vector<unsigned int>::const_iterator SingleCPUParticleData::end_types() const {
                    return cend_types();
                }

                std::vector<unsigned int>::const_iterator SingleCPUParticleData::cend_types() const {
                    return type->cbegin()+deactivated_index;
                }

                bool SingleCPUParticleData::isMarkedForDeactivation(const size_t index) {
                    return (*deactivated)[index];
                }

                size_t SingleCPUParticleData::getDeactivatedIndex() const {
                    return deactivated_index;
                }

                size_t SingleCPUParticleData::getNDeactivated() const {
                    return n_deactivated;
                }

                void SingleCPUParticleData::setParticleData(const readdy::model::Particle& particle, const size_t& index) {
                    (*ids)[index] = particle.getId();
                    (*positions)[index] = particle.getPos();
                    (*type)[index] = particle.getType();
                }

            }
        }
    }
}




