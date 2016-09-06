#include <readdy/model/Particle.h>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

/**
 * << detailed description >>
 *
 * @file Particle.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

namespace readdy {
    namespace model {
        Particle::Particle() : id(boost::uuids::random_generator()()), pos(0, 0, 0) {
        }

        bool Particle::operator==(const Particle &rhs) {
            return rhs.id == id;
        }

        bool Particle::operator!=(const Particle &rhs) {
            return !(*this == rhs);
        }

        Particle::Particle(double x, double y, double z, uint type) : id(boost::uuids::random_generator()()), pos(x, y, z) {
            this->type = type;
        }

        const Vec3 &Particle::getPos() const {
            return pos;
        }

        void Particle::setPos(const Vec3 &pos) {
            Particle::pos = pos;
        }

        unsigned int Particle::getType() const {
            return type;
        }

        void Particle::setType(unsigned int type) {
            Particle::type = type;
        }

        const boost::uuids::uuid &Particle::getId() const {
            return id;
        }

        Particle::Particle(Vec3 pos, unsigned int type, boost::uuids::uuid id) : pos(pos), type(type), id(id) {
        }

        void Particle::setId(const boost::uuids::uuid &id) {
            Particle::id = id;
        }

        Vec3 &Particle::getPos() {
            return pos;
        }


        Particle::~Particle() = default;

        std::ostream &operator<<(std::ostream &os, const Particle &p) {
            os << "Particle(id=" << boost::uuids::to_string(p.id) << ", type=" << p.type << ", pos=" << p.pos << ")";
            return os;
        }

    }
}
