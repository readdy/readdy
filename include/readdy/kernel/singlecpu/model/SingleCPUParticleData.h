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
#include <mutex>
#include <atomic>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace model {

class SingleCPUParticleData {
public:
    using marked_count_t = std::atomic<std::size_t>;

    // ctor / dtor
    SingleCPUParticleData();

    SingleCPUParticleData(bool useMarkedSet);

    SingleCPUParticleData(unsigned int capacity);

    SingleCPUParticleData(unsigned int capacity, bool useMarkedSet);

    ~SingleCPUParticleData();

    // move
    SingleCPUParticleData(SingleCPUParticleData &&rhs);

    SingleCPUParticleData &operator=(SingleCPUParticleData &&rhs);

    // copy
    SingleCPUParticleData(const SingleCPUParticleData &rhs) = delete;

    SingleCPUParticleData &operator=(const SingleCPUParticleData &rhs) = delete;

    std::vector<readdy::model::Particle::id_type>::iterator begin_ids();

    std::vector<readdy::model::Particle::id_type>::const_iterator begin_ids() const;

    std::vector<readdy::model::Particle::id_type>::const_iterator cbegin_ids() const;

    std::vector<readdy::model::Particle::id_type>::iterator end_ids();

    std::vector<readdy::model::Particle::id_type>::const_iterator end_ids() const;

    std::vector<readdy::model::Particle::id_type>::const_iterator cend_ids() const;

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

    std::vector<char>::iterator begin_deactivated();

    std::vector<char>::const_iterator begin_deactivated() const;

    std::vector<char>::const_iterator cbegin_deactivated() const;

    std::vector<char>::iterator end_deactivated();

    std::vector<char>::const_iterator end_deactivated() const;

    std::vector<char>::const_iterator cend_deactivated() const;

    void swap(SingleCPUParticleData &rhs);

    size_t size() const;

    size_t max_size() const;

    bool empty() const;

    void clear();

    void addParticle(const readdy::model::Particle &particle);

    void addParticles(const std::vector<readdy::model::Particle> &particles);

    /**
     * Remove a particle via its unique id.
     * @param particle the particle to be removed
     */
    void removeParticle(const readdy::model::Particle &particle);

    void removeParticle(const size_t index);

    void setParticleData(const readdy::model::Particle &particle, const size_t &index);

    readdy::model::Particle operator[](const size_t index) const;

    bool isMarkedForDeactivation(const size_t index);

    size_t getDeactivatedIndex() const;

    size_t getNDeactivated() const;

    void markForDeactivation(size_t index);

    template<typename T>
    void updateDeactivated(const T &update) const {
        std::lock_guard<std::mutex> lock(markedForDeactivationMutex);
        markedForDeactivation->insert(std::begin(update), std::end(update));
    }

    /**
     * This method is the counterpart to markForDeactivation.
     * The particles that were marked are now deactivated, i.e.,
     * for each marked particle:
     *   - If it is at the very end of the particle list, the
     *     counters are updated.
     *   - If not, the particle is swapped with the last active particle,
     *     so that again, all deactivated particles reside at the end
     *     of the internal data structure.
     */
    void deactivateMarked();

protected:
    std::unique_ptr<std::vector<readdy::model::Particle::id_type>> ids;
    std::unique_ptr<std::vector<readdy::model::Vec3>> positions;
    std::unique_ptr<std::vector<readdy::model::Vec3>> forces;
    std::unique_ptr<std::vector<unsigned int>> type;
    std::unique_ptr<std::vector<char>> deactivated;
    size_t deactivated_index;
    size_t n_deactivated;
    mutable marked_count_t n_marked;
    mutable std::unique_ptr<std::set<size_t>> markedForDeactivation;
    mutable std::mutex markedForDeactivationMutex;
    bool useMarkedSet = true;

    void deactivateMarkedSet();

    void deactivateMarkedNoSet();
};

}
}
}
}

#endif //READDY_MAIN_SINGLECPUPARTICLEDATA_H
