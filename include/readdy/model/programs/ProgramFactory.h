//
// Created by clonker on 08.04.16.
//

#ifndef READDY_MAIN_PROGRAMFACTORY_H
#define READDY_MAIN_PROGRAMFACTORY_H

#include <type_traits>
#include <unordered_map>
#include <vector>
#include <readdy/model/programs/Program.h>

namespace readdy {
    namespace model {
        namespace programs {
            class ProgramFactory {
            public:
                std::vector<std::string> getAvailablePrograms() const {
                    std::vector<std::string> names {};
                    for(auto&& e : factory) {
                        names.push_back(e.first);
                    }
                    return names;
                }

                template<typename T>
                std::unique_ptr<T> createProgramAs(const std::string &name) const {
                    auto&& it = factory.find(name);
                    if(it != factory.end()) {
                        return std::unique_ptr<T>(dynamic_cast<T *>(it->second()));
                    }
                    throw std::runtime_error("Could not find requested potential \""+name+"\" in factory.");
                }

                std::unique_ptr<Program> createProgram(std::string name) const {
                    return createProgramAs<Program>(name);
                }

            protected:
                std::unordered_map<std::string, std::function<Program *()>> factory;
            };
        }
    }
}

#endif //READDY_MAIN_PROGRAMFACTORY_H
