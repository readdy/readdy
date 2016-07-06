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
#include <boost/uuid/uuid.hpp>
#include "Vec3.h"

namespace readdy {
    namespace model {
        class Particle {
        public:
            const Vec3 &getPos() const;

            Vec3 &getPos();

            void setPos(const Vec3 &pos);

            unsigned int getType() const;

            void setType(const unsigned int type);

            const boost::uuids::uuid &getId() const;

            void setId(const boost::uuids::uuid &id);

            Particle();

            Particle(double x, double y, double z, unsigned int type);

            Particle(Vec3 pos, unsigned int type, boost::uuids::uuid id);

            virtual ~Particle();

            bool operator==(const Particle &rhs);

            bool operator!=(const Particle &rhs);


        private:
            Vec3 pos;
            unsigned int type;
            boost::uuids::uuid id;
        };
    }
}
#endif //READDY_MAIN_PARTICLE_H
