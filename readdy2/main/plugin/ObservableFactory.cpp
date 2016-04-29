#include <readdy/plugin/_internal/ObservableFactory.h>
#include <readdy/plugin/Observables.h>

readdy::plugin::_internal::ObservableFactory::ObservableFactory() {
    factory["ParticlePosition"] = [] { return new ParticlePositionObservable(); };
}

void readdy::plugin::_internal::ObservableFactory::registerObservable(const std::string &name, const std::function<Observable*()> create) {
    factory[name] = create;
}

std::unique_ptr<readdy::plugin::Observable> readdy::plugin::_internal::ObservableFactory::create(const std::string &name) {
    return std::unique_ptr<readdy::plugin::Observable>(factory[name]());
}

std::vector<std::string> readdy::plugin::_internal::ObservableFactory::getRegisteredObservableNames() const {
    std::vector<std::string> result;
    for(auto&& item : factory) {
        result.push_back(item.first);
    }
    return result;
}

