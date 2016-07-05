/**
 * << detailed description >>
 *
 * @file SimulationImpl.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.07.16
 */

struct Simulation::Impl {
    std::unique_ptr<readdy::model::Kernel> kernel;
    std::vector<std::unique_ptr<readdy::model::ObservableBase>> foo {};
    std::unordered_map<boost::uuids::uuid, std::unique_ptr<readdy::model::ObservableBase>, boost::hash<boost::uuids::uuid>> observables {};
    std::unordered_map<boost::uuids::uuid, boost::signals2::scoped_connection, boost::hash<boost::uuids::uuid>> observableConnections {};
};


template<typename T, typename... Args>
boost::uuids::uuid Simulation::registerObservable(const std::function<void(typename T::result_t)> callbackFun, unsigned int stride, Args... args) {
    ensureKernelSelected();
    boost::uuids::random_generator uuid_gen;
    auto uuid = uuid_gen();
    auto && obs = pimpl->kernel->createObservable<T>(stride, std::forward<Args>(args)...);
    obs->setCallback(std::move(callbackFun));
    pimpl->observables.emplace(uuid, std::move(obs));
    auto&& connection = pimpl->kernel->connectObservable(pimpl->observables[uuid].get());
    pimpl->observableConnections.emplace(uuid, std::move(connection));
    return uuid;
}
