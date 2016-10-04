/**
 * This header file contains the declaration of fission reactions, i.e., A->B+C. The products can be placed at a
 * specified distance and weighted by two members.
 *
 * @file Fission.h
 * @brief Declaration of fission reactions, i.e., A->B+C.
 * @author clonker
 * @date 20.06.16
 */

#ifndef READDY_MAIN_FISSION_H
#define READDY_MAIN_FISSION_H

#include "Reaction.h"
#include <boost/log/trivial.hpp>

namespace readdy {
namespace model {
namespace reactions {

class Fission : public Reaction<1> {

public:
    Fission(const std::string &name, unsigned int from, unsigned int to1, unsigned int to2,
            const double rate, const double productDistance, const double weight1 = 0.5,
            const double weight2 = 0.5) :
            Reaction(name, rate, 0, productDistance, 2), weight1(weight1), weight2(weight2) {
        educts = {from};
        products = {to1, to2};
        const auto sum = weight1 + weight2;
        if (sum != 1) {
            this->weight1 /= sum;
            this->weight2 /= sum;
            BOOST_LOG_TRIVIAL(warning) <<
                                       "The weights did not add up to 1, they were changed to weight1=" <<
                                       this->weight1 << ", weight2=" << this->weight2;
        }
    }

    const unsigned int getFrom() const {
        return educts[0];
    }

    const unsigned int getTo1() const {
        return products[0];
    }

    const unsigned int getTo2() const {
        return products[1];
    }

    const double getWeight1() const {
        return weight1;
    }

    const double getWeight2() const {
        return weight2;
    }

    virtual Fission *replicate() const override {
        return new Fission(*this);
    }

protected:
    double weight1, weight2;


};
}
}
}
#endif //READDY_MAIN_FISSION_H
