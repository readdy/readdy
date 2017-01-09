/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


#include <readdy/model/Particle.h>
#include <ostream>

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

std::atomic<Particle::id_type> Particle::id_counter {0};

Particle::Particle() : Particle(0,0,0,0) {
}

bool Particle::operator==(const Particle &rhs) const{
    return rhs.id == id;
}

bool Particle::operator!=(const Particle &rhs) const{
    return !(*this == rhs);
}

Particle::Particle(double x, double y, double z, unsigned int type) : id(std::atomic_fetch_add<unsigned long>(&id_counter, 1L)), pos(x, y, z) {
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

const Particle::id_type Particle::getId() const {
    return id;
}

Particle::Particle(Vec3 pos, unsigned int type, id_type id) : pos(pos), type(type), id(id) {
}

void Particle::setId(const id_type id) {
    Particle::id = id;
}

Vec3 &Particle::getPos() {
    return pos;
}

Particle::Particle(Vec3 pos, unsigned int type) : pos(pos), type(type), id(std::atomic_fetch_add<id_type>(&id_counter, 1)) {
}


Particle::~Particle() = default;

std::ostream &operator<<(std::ostream &os, const Particle &p) {
    os << "Particle(id=" << p.id << ", type=" << p.type << ", pos=" << p.pos << ")";
    return os;
}

}
}
