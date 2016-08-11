/**
 * This header file contains the definition of the PotentialFactory. It contains a map of available potentials to
 * a std function object that will create a new instance upon invocation.
 *
 * @file PotentialFactory.h
 * @brief Header file containing the definition of the PotentialFactory.
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIALFACTORY_H
#define READDY_MAIN_POTENTIALFACTORY_H

#include <unordered_map>
#include <functional>

#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/potentials/PotentialsOrder2.h>

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
                std::unique_ptr<T> createPotentialAs(const std::string &name) const {
                    auto&& it = factory.find(name);
                    if(it != factory.end()) {
                        return std::unique_ptr<T>(dynamic_cast<T *>(it->second()));
                    }
                    throw std::runtime_error("Could not find requested potential \""+name+"\" in factory.");
                }

                template<typename T>
                std::unique_ptr<T> createPotential() const {
                    return createPotentialAs<T>(getPotentialName<T>());
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
