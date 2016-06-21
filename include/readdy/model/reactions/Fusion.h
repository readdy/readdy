/**
 * << detailed description >>
 *
 * @file Fusion.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#ifndef READDY_MAIN_FUSION_H
#define READDY_MAIN_FUSION_H

#include "Reaction.h"

namespace readdy {
    namespace model {
        namespace reactions {
            class Fusion : public Reaction<2> {

            public:
                Fusion(const std::string &name, unsigned int from1, unsigned int from2, unsigned int to, const double &rate, const double &eductDistance) :
                        Reaction(name, rate, eductDistance, 0, 1)
                {
                    educts = {from1, from2};
                    products = {to};
                }

                const double getEductDistance() const {
                    return eductDistance;
                }

                const unsigned int getFrom1() const {
                    return educts[0];
                }

                const unsigned int getFrom2() const {
                    return educts[1];
                }

                const unsigned int getTo() const {
                    return products[0];
                }

            };
        }
    }
}

#endif //READDY_MAIN_FUSION_H
