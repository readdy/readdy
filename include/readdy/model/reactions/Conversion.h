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
            class Conversion : public Reaction {

            public:
                Conversion(const std::string &name, unsigned int typeFrom, unsigned int typeTo, const double &rate) :
                        Reaction(name, rate), typeFrom(typeFrom), typeTo(typeTo)
                { }

                unsigned int getTypeFrom() const {
                    return typeFrom;
                }

                unsigned int getTypeTo() const {
                    return typeTo;
                }

            protected:
                const unsigned int typeFrom, typeTo;
            };
        }
    }
}
#endif //READDY_MAIN_CONVERSION_H
