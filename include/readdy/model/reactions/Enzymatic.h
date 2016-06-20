/**
 * << detailed description >>
 *
 * @file Enzymatic.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#ifndef READDY_MAIN_ENZYMATIC_H
#define READDY_MAIN_ENZYMATIC_H

#include "Reaction.h"

namespace readdy {
    namespace model {
        namespace reactions {
            class Enzymatic : public Reaction {

            public:
                Enzymatic(const std::string &name, unsigned int catalyst, unsigned int from, unsigned int to, const double &rate, const double &radius) :
                        Reaction(name, rate), catalyst(catalyst), from(from), to(to), radius(radius) { }


                unsigned int getCatalyst() const {
                    return catalyst;
                }

                unsigned int getFrom() const {
                    return from;
                }

                unsigned int getTo() const {
                    return to;
                }

                double getRadius() const {
                    return radius;
                }

            protected:
                const unsigned int catalyst, from, to;
                const double radius;
            };
        }
    }
}
#endif //READDY_MAIN_ENZYMATIC_H
