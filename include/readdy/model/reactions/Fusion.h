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
            class Fusion : public Reaction {

            public:
                Fusion(const std::string &name, unsigned int from1, unsigned int from2, unsigned int to, const double &rate, const double &eductDistance) :
                        Reaction(name, rate), from1(from1), from2(from2), to(to), eductDistance(eductDistance) { }

                const double getEductDistance() const {
                    return eductDistance;
                }

                unsigned int getFrom1() const {
                    return from1;
                }

                unsigned int getFrom2() const {
                    return from2;
                }

                unsigned int getTo() const {
                    return to;
                }

            protected:
                const unsigned int from1, from2, to;
                const double eductDistance;
            };
        }
    }
}

#endif //READDY_MAIN_FUSION_H
