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

#ifndef READDY_MAIN_OBSERVABLEFACTORY_H
#define READDY_MAIN_OBSERVABLEFACTORY_H

#include <string>
#include <unordered_map>
#include <readdy/model/Observable.h>
#include <readdy/common/Utils.h>
#include <readdy/model/Observables.h>

namespace readdy {
    namespace model {
        class Kernel;
        namespace _internal {
            class NoSuchObservableException : public std::runtime_error {
            public:
                NoSuchObservableException(const std::string &__arg);
            };

            class ObservableFactory {
            public:
                ObservableFactory(Kernel *const kernel);

                void registerObservable(const std::string &name, const std::function<readdy::model::ObservableBase *()> create);

                std::unique_ptr<ObservableBase> create(const std::string &name);

                template<typename T, typename Obs1, typename Obs2>
                inline std::unique_ptr<T> create(Obs1 *obs1, Obs2 *obs2, unsigned int stride = 1) {
                    return std::make_unique<T>(kernel, obs1, obs2, stride);
                };

                template<typename T>
                inline std::unique_ptr<T> create() {
                    const auto name = ObservableName<T>::value;
                    if (readdy::utils::collections::hasKey(factory, name)) {
                        return std::unique_ptr<T>(dynamic_cast<T *>(factory[name]()));
                    }
                    throw NoSuchObservableException("The requested observable \"" + std::string(name) + "\" was not registered in the observable factory.");
                }

                std::vector<std::string> getRegisteredObservableNames() const;

            private:
                std::unordered_map<std::string, std::function<ObservableBase *()>> factory;
                Kernel *const kernel;
            };
        }
    }
}
#endif //READDY_MAIN_OBSERVABLEFACTORY_H
