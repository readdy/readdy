/**
 * << detailed description >>
 *
 * @file Fission.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#ifndef READDY_MAIN_FISSION_H
#define READDY_MAIN_FISSION_H

#include "Reaction.h"

namespace readdy {
    namespace model {
        namespace reactions {
            class Fission : public Reaction<1> {

            public:
                Fission(const std::string &name, unsigned int from, unsigned int to1, unsigned int to2, const double productDistance, const double &rate) :
                        Reaction(name, rate, 0, productDistance, 2)
                {
                    educts = {from};
                    products = {to1, to2};
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

                const double getProductDistance() const {
                    return productDistance;
                }

            };
        }
    }
}
#endif //READDY_MAIN_FISSION_H
