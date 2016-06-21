/**
 * << detailed description >>
 *
 * @file Death.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */


#ifndef READDY_MAIN_DEATH_H
#define READDY_MAIN_DEATH_H

#include "Reaction.h"

namespace readdy {
    namespace model {
        namespace reactions {
            class Death : public Reaction<1> {

            public:
                Death(const std::string &name, unsigned int typeFrom, const double &rate) : Reaction(name, rate, 0, 0, 0) {
                    educts[0] = typeFrom;
                }

                const unsigned int getTypeFrom() const {
                    return educts[0];
                }
            };
        }
    }
}
#endif //READDY_MAIN_DEATH_H
