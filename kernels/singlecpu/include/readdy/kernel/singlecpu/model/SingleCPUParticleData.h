/**
 * << detailed description >>
 *
 * @file SingleCPUParticleData.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#ifndef READDY_MAIN_SINGLECPUPARTICLEDATA_H
#define READDY_MAIN_SINGLECPUPARTICLEDATA_H

#include <memory>
#include <vector>
#include <readdy/model/Particle.h>
#include <set>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace model {

                class SingleCPUParticleData {
                public:

                    // ctor / dtor
                    SingleCPUParticleData();

                    SingleCPUParticleData(unsigned int capacity);

                    ~SingleCPUParticleData();

                    // move
                    SingleCPUParticleData(SingleCPUParticleData &&rhs);

                    SingleCPUParticleData &operator=(SingleCPUParticleData &&rhs);

                    // copy
                    SingleCPUParticleData(const SingleCPUParticleData &rhs) = delete;

                    SingleCPUParticleData &operator=(const SingleCPUParticleData &rhs) = delete;

                     std::vector<boost::uuids::uuid>::iterator begin_ids();
                     std::vector<boost::uuids::uuid>::const_iterator begin_ids() const;
                     std::vector<boost::uuids::uuid>::const_iterator cbegin_ids() const;
                     std::vector<boost::uuids::uuid>::iterator end_ids();
                     std::vector<boost::uuids::uuid>::const_iterator end_ids() const;
                     std::vector<boost::uuids::uuid>::const_iterator cend_ids() const;

                     std::vector<readdy::model::Vec3>::iterator begin_positions();
                     std::vector<readdy::model::Vec3>::const_iterator begin_positions() const;
                     std::vector<readdy::model::Vec3>::const_iterator cbegin_positions() const;
                     std::vector<readdy::model::Vec3>::iterator end_positions();
                     std::vector<readdy::model::Vec3>::const_iterator end_positions() const;
                     std::vector<readdy::model::Vec3>::const_iterator cend_positions() const;

                     std::vector<readdy::model::Vec3>::iterator begin_forces();
                     std::vector<readdy::model::Vec3>::const_iterator begin_forces() const;
                     std::vector<readdy::model::Vec3>::const_iterator cbegin_forces() const;
                     std::vector<readdy::model::Vec3>::iterator end_forces();
                     std::vector<readdy::model::Vec3>::const_iterator end_forces() const;
                     std::vector<readdy::model::Vec3>::const_iterator cend_forces() const;

                     std::vector<unsigned int>::iterator begin_types();
                     std::vector<unsigned int>::const_iterator begin_types() const;
                     std::vector<unsigned int>::const_iterator cbegin_types() const;
                     std::vector<unsigned int>::iterator end_types();
                     std::vector<unsigned int>::const_iterator end_types() const;
                     std::vector<unsigned int>::const_iterator cend_types() const;

                    void swap(SingleCPUParticleData &rhs);

                    size_t size();

                    size_t max_size();

                    bool empty();

                    void addParticle(const readdy::model::Particle &particle);
                    void addParticles(const std::vector<readdy::model::Particle> &particles);

                    void removeParticle(const readdy::model::Particle &particle);

                    void removeParticle(const size_t index);

                    readdy::model::Particle operator[](const size_t index);

                    bool isMarkedForDeactivation(int index);
                    size_t getDeactivatedIndex() const;
                    size_t getNDeactivated() const;

                    void markForDeactivation(size_t index);
                    void deactivateMarked();

                protected:
                    std::unique_ptr<std::vector<boost::uuids::uuid>> ids;
                    std::unique_ptr<std::vector<readdy::model::Vec3>> positions;
                    std::unique_ptr<std::vector<readdy::model::Vec3>> forces;
                    std::unique_ptr<std::vector<unsigned int>> type;
                    std::unique_ptr<std::vector<bool>> deactivated;
                    std::unique_ptr<std::set<size_t>> markedForDeactivation;
                    size_t deactivated_index;
                    size_t n_deactivated;
                };

            }
        }
    }
}

#endif //READDY_MAIN_SINGLECPUPARTICLEDATA_H
