/**
 * Decay reactions simply remove particles of a certain type with a certain rate of the system.
 *
 * @file Death.h
 * @brief This file contains the declaration of decay reactions.
 * @author clonker
 * @date 21.06.16
 */


#ifndef READDY_MAIN_DEATH_H
#define READDY_MAIN_DEATH_H

#include "Reaction.h"

namespace readdy {
namespace model {
namespace reactions {

class Decay : public Reaction<1> {

public:
    Decay(const std::string &name, unsigned int typeFrom, const double rate) : Reaction(name, rate, 0, 0, 0) {
        educts[0] = typeFrom;
    }

    const unsigned int getTypeFrom() const {
        return educts[0];
    }

    virtual void perform(const Particle &p1_in, const Particle &p2_in, Particle &p1_out, Particle &p2_out,
                         const rnd_ptr &rnd = nullptr) const override {};

    virtual Decay *replicate() const override {
        return new Decay(*this);
    }
};
}
}
}
#endif //READDY_MAIN_DEATH_H
