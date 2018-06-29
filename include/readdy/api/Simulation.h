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
 * This file contains mainly the simulation class. It can be
 * instantiated by using the default constructor and provides means to start a
 * simulation in a relatively simple manner. However, usage of the python top level API is recommended.
 *
 * @file Simulation.h
 * @brief File containing the ReaDDy top level API.
 * @author clonker
 * @date 08.02.16
 */

#pragma once

#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/Kernel.h>
#include <readdy/api/SimulationLoop.h>
#include <readdy/model/topologies/reactions/StructuralTopologyReaction.h>
#include "ObservableHandle.h"

NAMESPACE_BEGIN(readdy)
/**
 * Simulation is the focus of the high-level C++ API of ReaDDy.
 * This is where the system is set up and run for a certain number of
 * time steps.
 * Things like temperature, boxsize, reactions and potentials belong to the
 * context and are given to the kernel when run() is called.
 */
class Simulation {
public:
    /**
     * Type of an observable callback corresponding to a certain observable
     */
    template<typename T>
    using observable_callback = typename std::function<void(typename T::result_type)>;

    explicit Simulation(plugin::KernelProvider::kernel_ptr kernel) : _kernel(std::move(kernel)) {}

    /**
     * The default constructor.
     */
    explicit Simulation(const std::string &kernel) : Simulation(
            plugin::KernelProvider::getInstance().create(kernel)) {};

    /**
     * Creates a topology particle of a certain type at a position without adding it to the simulation box yet.
     * In order to instantiate it in the simulation it has to be used for creating a topology that owns this particular
     * particle instance.
     * @param type the type of topology flavored particle
     * @param pos the position
     * @return the created topology particle instance
     * @see addTopology()
     */
    readdy::model::TopologyParticle createTopologyParticle(const std::string &type, const Vec3 &pos) const {
        return _kernel->createTopologyParticle(type, pos);
    }

    /**
     * Checks whether the selected kernel supports topologies
     * @return true if the kernel supports topologies
     */
    bool kernelSupportsTopologies() const {
        return _kernel->supportsTopologies();
    }

    /**
     * Method for adding a new topology based on a topology type and topology particles. Node that initially the
     * particles are unconnected in the topology's connectivity graph, one has to make sure that there is a connection
     * through bonds between the particles upon simulation start.
     * @param type the type of topology to create
     * @param particles the topology's particles
     * @return a pointer to the instance of topology that was created
     * @see createTopologyParticle(type, pos)
     */
    readdy::model::top::GraphTopology *addTopology(const std::string &type,
                                                   const std::vector<readdy::model::TopologyParticle> &particles) {
        if (kernelSupportsTopologies()) {
            auto typeId = _kernel->context().topologyRegistry().idOf(type);
            return _kernel->stateModel().addTopology(typeId, particles);
        }
        throw std::logic_error("the selected kernel does not support topologies!");
    }

    /**
     * Method yielding a collection of pointers to the currently existing topologies in the simulation.
     * @return a vector of graph topology pointers
     */
    std::vector<readdy::model::top::GraphTopology *> currentTopologies() {
        return _kernel->stateModel().getTopologies();
    }

    /**
     * Method yielding a collection of particles that corresponds to the particles currently contained in the given
     * topology instance.
     * @param topology the topology instance
     * @return a vector of particles
     */
    std::vector<model::Particle> getParticlesForTopology(const model::top::GraphTopology &topology) const {
        return _kernel->stateModel().getParticlesForTopology(topology);
    }

    const model::observables::ObservableFactory &observe() const {
        return _kernel->observe();
    }

    /**
     * Registers a predefined observable with the kernel. A list of available observables can be obtained by
     * getAvailableObservables().
     * @tparam observable type
     * @param observable the observable
     * @return a uuid with which the observable is associated
     */
    template<typename T>
    ObservableHandle registerObservable(std::unique_ptr<T> observable, detail::is_observable_type<T> * = 0) {
        return registerObservable(std::move(observable), [](const typename T::result_type & /*unused*/) {});
    }

    /**
     * Registers a predefined observable with the kernel together with a callback.
     * A list of available observables can be obtained by getAvailableObservables().
     * @tparam T observable type
     * @param observable the observable instance
     * @param callback the callback
     * @return a observable handle that allows for post-hoc modification of the observable
     */
    template<typename T>
    ObservableHandle registerObservable(std::unique_ptr<T> observable, const observable_callback<T> &callback,
                                        detail::is_observable_type<T> * = 0) {
        auto connection = _kernel->connectObservable(observable.get());
        observable->callback() = callback;
        _observables.push_back(std::move(observable));
        _observableConnections.push_back(std::move(connection));
        return ObservableHandle{_observables.back().get()};
    }

    /**
     * A method to access the particle positions of a certain type.
     * @param type the type
     * @return a vector containing the particle positions
     */
    std::vector<Vec3> getParticlePositions(const std::string &type) {
        auto typeId = _kernel->context().particleTypes().idOf(type);
        const auto particles = _kernel->stateModel().getParticles();
        std::vector<Vec3> positions;
        positions.reserve(particles.size());
        for (const auto &p : particles) {
            if (p.getType() == typeId) {
                positions.push_back(p.getPos());
            }
        }
        return positions;
    }


    /**
     * A method to add a particle of a previously registered type to the system.
     * @param x the x coordinate
     * @param y the y coordinate
     * @param z the z coordinate
     * @param type the type of the particle
     */
    void addParticle(const std::string &type, scalar x, scalar y, scalar z) {
        const auto &s = context().boxSize();
        if (fabs(x) <= .5 * s[0] && fabs(y) <= .5 * s[1] && fabs(z) <= .5 * s[2]) {
            _kernel->stateModel().addParticle({x, y, z, context().particleTypes().idOf(type)});
        } else {
            log::error("particle position was not in bounds of the simulation box!");
        }
    }

    /**
     * Method that gives access to all the positions of all the particles in the system.
     * @return a vector containing all the positions
     */
    const std::vector<Vec3> getAllParticlePositions() const {
        return _kernel->stateModel().getParticlePositions();
    }

    /**
     * Yields a modifyable reference to the current context of this simulation. Can be used to set time-independent
     * properties.
     * @return the context
     */
    model::Context &context() {
        return _kernel->context();
    }

    /**
     * Yields a nonmodifyable reference to the current context of this simulation.
     * @return cref to the context
     */
    const model::Context &context() const {
        return _kernel->context();
    }

    /**
     * Method that returns the type of the currently selected kernel.
     * Raises if there is no kernel selected.
     * @return the kernel id string
     */
    const std::string &selectedKernelType() const {
        return _kernel->name();
    }

    /**
     * Runs the simulation using the default simulation loop for a given number of time steps and a given time step.
     * @param steps the number of steps
     * @param timeStep the time step
     * @see runScheme()
     */
    virtual void run(time_step_type steps, scalar timeStep) {
        createLoop(timeStep).run(steps);
    }

    api::SimulationLoop createLoop(scalar timeStep) {
        return api::SimulationLoop(_kernel.get(), timeStep, _performanceRoot);
    }

    /**
     * Checks if the kernel is running on single precision.
     * @return true if it is running on single precision
     */
    bool singlePrecision() const {
        return _kernel->singlePrecision();
    }

    /**
     * Checks if the kernel is running on double precision.
     * @return true if it is running on double precision
     */
    bool doublePrecision() const {
        return _kernel->doublePrecision();
    }

    /**
     * Access the root node of the performance measurement tree.
     * @return reference to root node of performance measurement
     */
    const util::PerformanceNode &performanceRoot() {
        return _performanceRoot;
    }

    /**
     * This method can be used to set kernel specific configuration. It takes a JSON string.
     * @param conf json
     */
    void setKernelConfiguration(const std::string &conf) {
        _kernel->context().setKernelConfiguration(conf);
    }

    model::StateModel &stateModel() {
        return _kernel->stateModel();
    }

    const model::StateModel &stateModel() const {
        return _kernel->stateModel();
    }

private:
    plugin::KernelProvider::kernel_ptr _kernel;
    std::vector<std::unique_ptr<readdy::model::observables::ObservableBase>> _observables{};
    std::vector<readdy::signals::scoped_connection> _observableConnections{};
    util::PerformanceNode _performanceRoot{"simulation", true};
};

NAMESPACE_END(readdy)
