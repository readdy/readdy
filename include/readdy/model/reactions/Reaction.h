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
#include <readdy/model/Particle.h>

namespace readdy {
    namespace model {
        namespace reactions {
            template<unsigned int N_EDUCTS>
            class Reaction {

            public:
                Reaction(const std::string &name, const double &rate, const double &eductDistance, const double &productDistance, const unsigned int n_products) :
                        name(name),
                        id(boost::uuids::random_generator()()),
                        rate(rate),
                        eductDistance(eductDistance),
                        productDistance(productDistance),
                        _n_products(n_products) { }

                virtual ~Reaction() = default;


                const std::string &getName() const {
                    return name;
                }

                const boost::uuids::uuid &getId() const {
                    return id;
                }

                const double &getRate() const {
                    return rate;
                }

                const unsigned int &getNEducts() const {
                    return _n_educts;
                }

                const unsigned int &getNProducts() const {
                    return _n_products;
                }

                const double &getEductDistance() const {
                    return eductDistance;
                }

                const double &getProductDistance() const {
                    return productDistance;
                }

                virtual void perform(const Particle &p1_in, const Particle &p2_in, Particle &p1_out, Particle &p2_out) const { };

                /*template<unsigned int _N>
                friend std::ostream &operator<<(std::ostream& os, const Reaction& reaction);*/

            protected:
                const unsigned int _n_educts = N_EDUCTS;
                const unsigned int _n_products;
                std::array<unsigned int, N_EDUCTS> educts;
                std::array<unsigned int, 2> products;
                const std::string name;
                const boost::uuids::uuid id;
                const double rate;
                const double eductDistance;
                const double productDistance;
            };

            /*template<unsigned int _N>
            inline std::ostream& operator<<(std::ostream& os, const Reaction<_N>& reaction) {
                os << "Reaction(\"" << reaction.name <<"\", N_Educts="<<_N<<", N_Products="<<reaction._n_products<<", (";
                for(int i = 0; i < _N; i++) {
                    if(i > 0) os << ",";
                    os << reaction.educts[i];
                }
                os <<") -> (";
                for(int i = 0; i < reaction._n_products; i++) {
                    if(i > 0) os << ",";
                    os << reaction.products[i];
                }
                os <<"), rate="<<reaction.rate<<", eductDist="<<reaction.eductDistance<<", prodDist="<<reaction.productDistance<<")";
                return os;
            }*/

        }
    }
}

#endif //READDY_MAIN_REACTION_H
