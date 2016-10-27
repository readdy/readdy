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
    using particle_type = readdy::model::Particle;
    using index_t = std::vector<particle_type::id_type>::size_type;

    struct Entry {
        particle_type::id_type id;
        particle_type::pos_type pos;
        readdy::model::Vec3 force;
        particle_type::type_type type;
        bool deactivated;
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

    void addParticles(const std::vector<particle_type> &particles);
/*
    void removeParticle(const particle_type &particle);

    void removeParticle(const index_t index);

    particle_type operator[](const index_t index) const;

    bool isMarkedForDeactivation(const index_t index);

    index_t getNDeactivated() const;

    void markForDeactivation(const size_t index);
*/
protected:
    std::stack<index_t> blanks;
};

}
}
}
}

#endif //READDY_KERNEL_CPU_PARTICLEDATA_H
