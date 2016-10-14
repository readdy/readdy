#include <readdy/model/reactions/Reaction.h>
#include <readdy/model/potentials/Potential.h>
#include <readdy/common/signals.h>

/**
 * This file contains the gory implementation details of the Simulation.h template code and thus is separated
 * from the rest.
 *
 * @file SimulationImpl.tpp
 * @brief Header-implementation file for Simulation.h
 * @author clonker
 * @date 05.07.16
 */

struct Simulation::Impl {
    std::unique_ptr<readdy::model::Kernel> kernel;
    std::vector<std::unique_ptr<readdy::model::ObservableBase>> foo {};
    std::unordered_map<unsigned long, std::unique_ptr<readdy::model::ObservableBase>> observables {};
    std::unordered_map<unsigned long, readdy::signals::scoped_connection> observableConnections {};
    unsigned long counter = 0;
};


template<typename T, typename... Args>
unsigned long Simulation::registerObservable(const std::function<void(typename T::result_t)> callbackFun, unsigned int stride, Args... args) {
    ensureKernelSelected();
    auto uuid = pimpl->counter++;
    auto && obs = pimpl->kernel->createObservable<T>(stride, std::forward<Args>(args)...);
    obs->setCallback(std::move(callbackFun));
    pimpl->observables.emplace(uuid, std::move(obs));
    auto&& connection = pimpl->kernel->connectObservable(pimpl->observables[uuid].get());
    pimpl->observableConnections.emplace(uuid, std::move(connection));
    return uuid;
}


template<typename SchemeType>
readdy::api::SchemeConfigurator<SchemeType> Simulation::runScheme(bool useDefaults) {
    ensureKernelSelected();
    return readdy::api::SchemeConfigurator<SchemeType>(getSelectedKernel(), useDefaults);
}