/**
 * << detailed description >>
 *
 * @file ParticleData.h
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#ifndef READDY_KERNEL_CPU_PARTICLEDATA_H
#define READDY_KERNEL_CPU_PARTICLEDATA_H

#include <cstddef>
#include <atomic>
#include <vector>
#include <memory>
#include <stack>

#include <readdy/model/Particle.h>
#include <readdy/model/KernelContext.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace model {

class ParticleData {
    friend class NeighborList;
public:

    struct Entry;
    using Neighbor = std::size_t;
    using ctx_t = readdy::model::KernelContext;
    using particle_type = readdy::model::Particle;
    using entries_t = std::vector<Entry>;
    using neighbors_t = std::vector<Neighbor>;
    using neighbor_list_t = std::vector<neighbors_t>;
    using index_t = entries_t::size_type;
    using update_t = std::pair<ParticleData::entries_t, std::vector<index_t>>;
    using force_t = readdy::model::Vec3;
    using displacement_t = double;

    /**
     * Particle data entry with padding such that it fits exactly into 64 bytes.
     */
    struct Entry {
        Entry(const particle_type &particle) : id(particle.getId()), pos(particle.getPos()), force(force_t()),
                                               type(particle.getType()), deactivated(false), displacement(0) { }

        Entry(const Entry&) = delete;
        Entry& operator=(const Entry&) = delete;
        Entry(Entry&&);
        Entry& operator=(Entry&&);

        bool is_deactivated() const;
        const particle_type::pos_type &position() const;

        particle_type::id_type id; // 8 bytes
        force_t force; // 3*8 + 8 = 32 bytes
        displacement_t displacement; // 32 + 8 = 40 bytes

    private:
        friend class readdy::kernel::cpu_dense::model::ParticleData;
        friend class NeighborList;

        particle_type::pos_type pos; // 40 + 3*8 = 64 bytes
    public:
        particle_type::type_type type; // 64 + 4 = 68 bytes
    private:
        bool deactivated; // 68 + 1 = 69 bytes
        char padding[3]; // 69 + 3 = 72 bytes
    };

    // ctor / dtor
    ParticleData(readdy::model::KernelContext *const context);

    ~ParticleData();

    // delete move and copy
    ParticleData(ParticleData &&rhs) = delete;

    ParticleData &operator=(ParticleData &&rhs) = delete;

    ParticleData(const ParticleData &rhs) = delete;

    ParticleData &operator=(const ParticleData &rhs) = delete;

    std::size_t size() const;

    bool empty() const;

    void clear();

    void addParticle(const particle_type &particle);

    index_t addEntry(Entry &&entry);

    void addParticles(const std::vector<particle_type> &particles);

    readdy::model::Particle getParticle(const index_t index) const;

    readdy::model::Particle toParticle(const Entry& e) const;

    void removeParticle(const particle_type &particle);

    void removeParticle(const index_t index);

    void removeEntry(index_t entry);

    auto begin() -> decltype(std::declval<entries_t>().begin()) {
        return entries.begin();
    }
    auto end() -> decltype(std::declval<entries_t>().end()) {
        return entries.end();
    }
    auto cbegin() const -> decltype(std::declval<entries_t>().cbegin()) {
        return entries.cbegin();
    }
    auto cend() const -> decltype(std::declval<entries_t>().cend()) {
        return entries.cend();
    }
    auto begin() const -> decltype(std::declval<entries_t>().cbegin()) {
        return cbegin();
    }
    auto end() const -> decltype(std::declval<entries_t>().cend()) {
        return cend();
    }

    Entry& entry_at(index_t);
    const Entry& entry_at(index_t) const;
    const Entry& centry_at(index_t) const;

    neighbors_t& neighbors_at(index_t);
    const neighbors_t& neighbors_at(index_t) const;
    const neighbors_t& cneighbors_at(index_t) const;

    const particle_type::pos_type& pos(index_t) const;

    void setFixPosFun(const ctx_t::fix_pos_fun&);

    index_t getEntryIndex(const Entry *const entry) const;

    index_t getNDeactivated() const;

    /**
     *
     * @return vector of new entries
     */
    std::vector<index_t> update(update_t&&);
    void displace(Entry&, const particle_type::pos_type& delta);

protected:


    std::stack<index_t, std::vector<index_t>> blanks;
    neighbor_list_t neighbors;
    entries_t entries;
    ctx_t::fix_pos_fun fixPos;
};

}
}
}
}

#endif //READDY_KERNEL_CPU_PARTICLEDATA_H
