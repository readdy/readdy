/**
 * << detailed description >>
 *
 * ... should not be instantiated directly
 *
 * @file ObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.04.16
 */

#ifndef READDY2_MAIN_OBSERVABLEFACTORY_H
#define READDY2_MAIN_OBSERVABLEFACTORY_H

#include <string>
#include <unordered_map>
#include <readdy/plugin/Observable.h>

namespace readdy {
    namespace plugin {
        namespace _internal {
            class ObservableFactory {
            public:
                ObservableFactory();
                void registerObservable(const std::string &name, const std::function<Observable*()> create);
                std::unique_ptr<Observable> create(const std::string &name);
                std::vector<std::string> getRegisteredObservableNames() const;
            private:
                std::unordered_map<std::string, std::function<Observable*()>> factory;
            };
        }
    }
}
#endif //READDY2_MAIN_OBSERVABLEFACTORY_H
