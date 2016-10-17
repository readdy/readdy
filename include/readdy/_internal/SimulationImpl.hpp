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
    std::unordered_map<unsigned long, readdy::signals::scoped_connection> observableConnections {};
    std::unordered_map<unsigned long, std::unique_ptr<readdy::model::ObservableBase>> observables {};
    std::unique_ptr<readdy::model::Kernel> kernel;
    unsigned long counter = 0;
};


template<typename T, typename... Args>
unsigned long Simulation::registerObservable(const std::function<void(typename T::result_t)>& callbackFun, unsigned int stride, Args... args) {
    ensureKernelSelected();
    auto uuid = pimpl->counter++;
    auto obs = pimpl->kernel->createObservable<T>(stride, std::forward<Args>(args)...);
    obs->setCallback(callbackFun);
    auto&& connection = pimpl->kernel->connectObservable(obs.get());
    pimpl->observables.emplace(uuid, std::move(obs));
    pimpl->observableConnections.emplace(uuid, std::move(connection));
    return uuid;
}


template<typename SchemeType>
readdy::api::SchemeConfigurator<SchemeType> Simulation::runScheme(bool useDefaults) {
    ensureKernelSelected();
    return readdy::api::SchemeConfigurator<SchemeType>(getSelectedKernel(), useDefaults);
}