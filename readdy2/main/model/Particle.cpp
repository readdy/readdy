#include <readdy/model/Particle.h>
#include <boost/uuid/uuid_generators.hpp>

/**
 * << detailed description >>
 *
 * @file Particle.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

readdy::model::Particle::Particle() : id(boost::uuids::random_generator()()), pos(0,0,0) {
}

bool readdy::model::Particle::operator==(const readdy::model::Particle &rhs) {
    return rhs.id == id;
}

bool readdy::model::Particle::operator!=(const readdy::model::Particle &rhs) {
    return !(*this == rhs);
}

readdy::model::Particle::Particle(double x, double y, double z, uint type) : id(boost::uuids::random_generator()()), pos(x,y,z){
    this->type = type;
}


readdy::model::Particle::~Particle() = default;

