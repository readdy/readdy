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

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

class ParticleData {
public:
    struct Entry;

    using particle_type = readdy::model::Particle;
    using entries_t = std::vector<Entry>;
    using index_t = entries_t::size_type;
    using update_t = std::pair<ParticleData::entries_t, std::vector<ParticleData::Entry *>>;

    /**
     * Particle data entry with padding such that it fits exactly into 64 bytes.
     */
    struct Entry {
        Entry(const particle_type &particle) : id(particle.getId()), pos(particle.getPos()),
                                               force(readdy::model::Vec3()), type(particle.getType()),
                                               deactivated(false) { }

        particle_type::id_type id; // 4 bytes
        particle_type::pos_type pos; // 4 + 24 = 28 bytes
        readdy::model::Vec3 force; // 28 + 24 = 52 bytes
        particle_type::type_type type; // 52 + 8 = 60 bytes

        bool is_deactivated() const {
            return deactivated;
        }

    private:
        friend class readdy::kernel::cpu::model::ParticleData;

        bool deactivated; // 60 + 1 = 61 bytes
        bool padding[3]; // 61 + 3 = 64 bytes
    };

    // ctor / dtor
    ParticleData();

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

    void addEntries(const std::vector<Entry> &entries);

    ParticleData::Entry* addEntry(Entry entry);

    void addParticles(const std::vector<particle_type> &particles);

    readdy::model::Particle getParticle(const index_t index) const;

    readdy::model::Particle toParticle(const Entry& e) const;

    void removeParticle(const particle_type &particle);

    void removeParticle(const index_t index);

    void removeEntry(Entry *const entry);

    index_t getEntryIndex(const Entry *const entry) const;

    index_t getNDeactivated() const;

    void update(update_t&&);

    entries_t entries;

protected:
    std::stack<index_t> blanks;
};

}
}
}
}

#endif //READDY_KERNEL_CPU_PARTICLEDATA_H
