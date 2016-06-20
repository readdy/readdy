/**
 * Base class for all possible types of reaction. Currently:
 *  - Conversion: A -> B
 *  - Enzymatic: A + C -> B + C where C is a catalyst
 *  - Fission: A -> B + C
 *  - Fusion: A + B -> C
 *
 * @file Reactions.h
 * @brief Reaction base class.
 * @author clonker
 * @date 17.06.16
 */

#ifndef READDY_MAIN_REACTION_H
#define READDY_MAIN_REACTION_H

#include <boost/uuid/uuid.hpp>
#include <string>
#include <boost/uuid/random_generator.hpp>

namespace readdy {
    namespace model {
        namespace reactions {
            class Reaction {

            public:
                Reaction(const std::string &name, const double &rate) :
                        name(name),
                        id(boost::uuids::random_generator()()),
                        rate(rate)
                { }
                virtual ~Reaction() = default;


                const std::string &getName() const {
                    return name;
                }

                const boost::uuids::uuid &getId() const {
                    return id;
                }

                const double getRate() const {
                    return rate;
                }

            protected:
                const std::string name;
                const boost::uuids::uuid id;
                const double rate;
            };
        }
    }
}

#endif //READDY_MAIN_REACTION_H
