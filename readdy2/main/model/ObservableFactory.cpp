#include <readdy/model/_internal/ObservableFactory.h>
#include <readdy/model/Observables.h>

namespace readdy {
    namespace model {
        namespace _internal {
            ObservableFactory::ObservableFactory(Kernel *const kernel) {
                factory["ParticlePosition"] = [kernel] { return new ParticlePositionObservable(kernel); };
            }

            void ObservableFactory::registerObservable(const std::string &name, const std::function<Observable *()> create) {
                factory[name] = create;
            }

            std::unique_ptr<Observable> ObservableFactory::create(const std::string &name) {
                return std::unique_ptr<Observable>(factory[name]());
            }

            std::vector<std::string> ObservableFactory::getRegisteredObservableNames() const {
                std::vector<std::string> result;
                for (auto &&item : factory) {
                    result.push_back(item.first);
                }
                return result;
            }


        }
    }
}
