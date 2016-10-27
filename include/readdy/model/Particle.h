/**
 * A particle is a composite of position, type, and id. The type is an "unsigned int" value, which is mapped by
 * the KernelContext.
 *
 * @file Particle.h
 * @brief Header file containing the definitions for the particle class.
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_PARTICLE_H
#define READDY_MAIN_PARTICLE_H

#include <array>
#include <string>
#include <atomic>
#include <spdlog/fmt/ostr.h>
#include "Vec3.h"

namespace readdy {
namespace model {

class Particle {
public:
    
    using id_type = unsigned long;
    using pos_type = Vec3;
    using type_type = unsigned int;
    
    const Vec3 &getPos() const;

    Vec3 &getPos();

    void setPos(const Vec3 &pos);

    unsigned int getType() const;

    void setType(const unsigned int type);

    const id_type getId() const;

    void setId(const id_type id);

    Particle();

    Particle(double x, double y, double z, unsigned int type);

    Particle(Vec3 pos, unsigned int type);
    Particle(Vec3 pos, unsigned int type, id_type id);

    virtual ~Particle();

    bool operator==(const Particle &rhs) const;

    bool operator!=(const Particle &rhs) const;

    friend std::ostream &operator<<(std::ostream &, const Particle &);


private:
    Vec3 pos;
    unsigned int type;
    id_type id;

    static std::atomic<id_type> id_counter;
};

}
}
#endif //READDY_MAIN_PARTICLE_H
