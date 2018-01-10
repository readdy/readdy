/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file Simulation_misc.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/plugin/KernelProvider.h>
#include "../Simulation.h"

NAMESPACE_BEGIN(readdy)

struct Simulation::Impl {
    plugin::KernelProvider::kernel_ptr kernel;
    std::unordered_map<unsigned long, std::unique_ptr<readdy::model::observables::ObservableBase>> observables{};
    std::unordered_map<unsigned long, readdy::signals::scoped_connection> observableConnections{};
    unsigned long counter = 0;
    util::PerformanceNode performanceRoot{"simulation", true};
};

template<typename T>
inline ObservableHandle Simulation::registerObservable(std::unique_ptr<T> observable, detail::is_observable_type<T>*) {
    ensureKernelSelected();
    auto uuid = pimpl->counter++;
    auto connection = pimpl->kernel->connectObservable(observable.get());
    pimpl->observables.emplace(uuid, std::move(observable));
    pimpl->observableConnections.emplace(uuid, std::move(connection));
    return {uuid, pimpl->observables.at(uuid).get()};
}

template<typename T>
inline ObservableHandle Simulation::registerObservable(std::unique_ptr<T> observable,
                                                       const observable_callback<T> &callback,
                                                       detail::is_observable_type<T>*) {
    ensureKernelSelected();
    auto uuid = pimpl->counter++;
    auto connection = pimpl->kernel->connectObservable(observable.get());
    observable->setCallback(callback);
    pimpl->observables.emplace(uuid, std::move(observable));
    pimpl->observableConnections.emplace(uuid, std::move(connection));
    return {uuid, pimpl->observables.at(uuid).get()};
}

template<typename SchemeType>
inline readdy::api::SchemeConfigurator<SchemeType> Simulation::runScheme(bool useDefaults) {
    ensureKernelSelected();
    return readdy::api::SchemeConfigurator<SchemeType>(getSelectedKernel(), pimpl->performanceRoot, useDefaults);
}
NAMESPACE_END(readdy)
