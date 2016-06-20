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
            class Fission : public Reaction {

            public:
                Fission(const std::string &name, unsigned int from, unsigned int to1, unsigned int to2, const double productDistance, const double &rate) :
                        Reaction(name, rate), from(from), to1(to1), to2(to2), productDistance(productDistance)
                { }

                unsigned int getFrom() const {
                    return from;
                }

                unsigned int getTo1() const {
                    return to1;
                }

                unsigned int getTo2() const {
                    return to2;
                }

                const double getProductDistance() const {
                    return productDistance;
                }

            protected:
                const unsigned int from, to1, to2;
                const double productDistance;
            };
        }
    }
}
#endif //READDY_MAIN_FISSION_H
