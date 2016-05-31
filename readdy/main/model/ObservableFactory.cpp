#include <readdy/model/_internal/ObservableFactory.h>
#include <readdy/model/Observables.h>

namespace readdy {
    namespace model {
        namespace _internal {
            ObservableFactory::ObservableFactory(Kernel *const kernel) : kernel(kernel) {
                factory[ObservableName<ParticlePositionObservable>::value] = [kernel] { return new ParticlePositionObservable(kernel); };
            }

            void ObservableFactory::registerObservable(const std::string &name, const std::function<ObservableBase *()> create) {
                factory[name] = create;
            }

            std::unique_ptr<ObservableBase> ObservableFactory::create(const std::string &name) {
                if (readdy::utils::collections::hasKey(factory, name)) {
                    return std::unique_ptr<ObservableBase>(factory[name]());
                }
                throw NoSuchObservableException("The requested observable \"" + name + "\" was not registered in the observable factory.");
            }

            std::vector<std::string> ObservableFactory::getRegisteredObservableNames() const {
                std::vector<std::string> result;
                for (auto &&item : factory) {
                    result.push_back(item.first);
                }
                return result;
            }

            NoSuchObservableException::NoSuchObservableException(const std::string &__arg) : std::runtime_error(__arg) {
            }

        }
    }
}


