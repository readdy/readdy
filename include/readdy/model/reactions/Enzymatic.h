/**
 * This file contains the declaration for enzymatic reactions, i.e., A+C->B+C.
 *
 * @file Enzymatic.h
 * @brief Declaration file of enzymatic reactions, i.e., A+C->B+C.
 * @author clonker
 * @date 20.06.16
 */

#ifndef READDY_MAIN_ENZYMATIC_H
#define READDY_MAIN_ENZYMATIC_H

#include "Reaction.h"

namespace readdy {
namespace model {
namespace reactions {

class Enzymatic : public Reaction<2> {

public:
    Enzymatic(const std::string &name, unsigned int catalyst, unsigned int from, unsigned int to, const double rate,
              const double eductDistance) :
            Reaction(name, rate, eductDistance, 0, 2) {
        educts = {from, catalyst};
        products = {to, catalyst};
    }


    const unsigned int getCatalyst() const {
        return educts[1];
    }

    const unsigned int getFrom() const {
        return educts[0];
    }

    const unsigned int getTo() const {
        return products[0];
    }

    virtual const ReactionType getType() override {
        return ReactionType::Enzymatic;
    }
};
}
}
}
#endif //READDY_MAIN_ENZYMATIC_H
