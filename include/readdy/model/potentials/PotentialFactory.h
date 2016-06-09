/**
 * << detailed description >>
 *
 * @file PotentialFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIALFACTORY_H
#define READDY_MAIN_POTENTIALFACTORY_H

#include <unordered_map>
#include <functional>

#include <readdy/model/potentials/Potential.h>

namespace readdy {
    namespace model {

        class Kernel;
        namespace potentials {
            class PotentialFactory {
            public:
                std::vector<std::string> getAvailablePotentials() const {
                    std::vector<std::string> names {factory.size()};
                    for(auto&& e : factory) {
                        names.push_back(e.first);
                    }
                    return names;
                }

                template<typename T>
                std::unique_ptr<T> createPotentialAs(std::string name) const {
                    return std::unique_ptr<T>(dynamic_cast<T*>(factory.find(name)->second()));
                }

                std::unique_ptr<Potential> createPotential(std::string name) const {
                    return createPotentialAs<Potential>(name);
                }

            protected:
                std::unordered_map<std::string, std::function<Potential *()>> factory;
            };

        }
    }
}

#endif //READDY_MAIN_POTENTIALFACTORY_H
