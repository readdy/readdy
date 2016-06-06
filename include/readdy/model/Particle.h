/**
 * << detailed description >>
 *
 * @file Particle.h
 * @brief << brief description >>
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

            void setPos(const Vec3 &pos);

            unsigned int getType() const;

            void setType(const unsigned int type);

            const boost::uuids::uuid &getId() const;

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
