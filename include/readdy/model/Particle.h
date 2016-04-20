/**
 * << detailed description >>
 *
 * @file Particle.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY2_MAIN_PARTICLE_H
#define READDY2_MAIN_PARTICLE_H

#include <array>
#include <string>
#include <boost/uuid/uuid.hpp>

namespace readdy {
    namespace model {
        class Particle {
        public:
            std::array<double, 3> pos;
            std::string type;

            Particle();
            Particle(double x, double y, double z, std::string type);
            virtual ~Particle();

            bool operator==(const Particle& rhs);
            bool operator!=(const Particle& rhs);

        private:
            boost::uuids::uuid id;
        };
    }
}
#endif //READDY2_MAIN_PARTICLE_H
