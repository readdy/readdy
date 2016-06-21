/**
 * << detailed description >>
 *
 * @file Conversion.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#ifndef READDY_MAIN_CONVERSION_H
#define READDY_MAIN_CONVERSION_H

#include "Reaction.h"

namespace readdy {
    namespace model {
        namespace reactions {
            class Conversion : public Reaction<1> {

            public:
                Conversion(const std::string &name, unsigned int typeFrom, unsigned int typeTo, const double &rate) :
                        Reaction(name, rate, 0 ,0, 1)
                {
                    educts = {typeFrom};
                    products = {typeTo};
                }

                const unsigned int getTypeFrom() const {
                    return educts[0];
                }

                const unsigned int getTypeTo() const {
                    return products[0];
                }

            };
        }
    }
}
#endif //READDY_MAIN_CONVERSION_H
