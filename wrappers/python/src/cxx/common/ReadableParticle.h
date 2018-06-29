//
// Created by mho on 5/29/18.
//

#pragma once

#include <utility>

#include <readdy/model/Particle.h>
#include <readdy/model/Context.h>

#include "common.h"

namespace rpy {

class ReadableParticle {
public:
    ReadableParticle(const readdy::model::Particle &particle, const readdy::model::Context &context)
            : _pos(3, particle.getPos().data.data()), _type(context.particleTypes().nameOf(particle.getType())),
              _id(particle.getId()) {}

    const auto &pos() const {
        return _pos;
    }

    const auto &type() const {
        return _type;
    }

    const auto &id() const {
        return _id;
    }

private:
    np_array<readdy::scalar> _pos;
    std::string _type;
    readdy::model::Particle::id_type _id;
};

}
