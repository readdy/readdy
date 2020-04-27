/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
#include <readdy/api/ObservableHandle.h>

namespace readdy {
/**
 * Simulation is the focus of the high-level C++ API of ReaDDy.
 * This is where the system is set up and run for a certain number of
 * time steps.
 * Things like temperature, boxsize, reactions and potentials belong to the
 * Context.
 */
class Simulation {
public:
    /**
     * Type of an observable callback corresponding to a certain observable
     */
    template<typename T>
    using observable_callback = typename std::function<void(typename T::result_type)>;

    /** User provided context is copied, and no modifiable view on it is provided after construction of Simulation */
    explicit Simulation(plugin::KernelProvider::kernel_ptr kernel, model::Context ctx) : _kernel(std::move(kernel)) {
        _kernel->context() = std::move(ctx);
    }

    /**
     * The default constructor.
     */
    explicit Simulation(const std::string &kernel, model::Context ctx) : Simulation(
            plugin::KernelProvider::getInstance().create(kernel), std::move(ctx)) {};

    /** Wrap a kernel with a working context in a simulation */
    explicit Simulation(plugin::KernelProvider::kernel_ptr kernel) : _kernel(std::move(kernel)) {}

    /**
     * Creates a topology particle of a certain type at a position without adding it to the simulation box yet.
     * In order to instantiate it in the simulation it has to be used for creating a topology that owns this particular
     * particle instance.
     * @param type the type of topology flavored particle
     * @param pos the position
     * @return the created topology particle instance
     * @see addTopology()
     */
    [[nodiscard]] readdy::model::Particle createTopologyParticle(const std::string &type, const Vec3 &pos) const {
        return _kernel->createTopologyParticle(type, pos);
    }

    /**
     * Checks whether the selected kernel supports topologies
     * @return true if the kernel supports topologies
     */
    [[nodiscard]] bool kernelSupportsTopologies() const {
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
                                                   const std::vector<readdy::model::Particle> &particles) {
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

    [[nodiscard]] auto currentParticles() const {
        return _kernel->stateModel().getParticles();
    }

    /**
     * Method yielding a collection of particles that corresponds to the particles currently contained in the given
     * topology instance.
     * @param topology the topology instance
     * @return a vector of particles
     */
    [[nodiscard]] std::vector<model::Particle> getParticlesForTopology(const model::top::GraphTopology &topology) const {
        return _kernel->stateModel().getParticlesForTopology(topology);
    }

    [[nodiscard]] const model::observables::ObservableFactory &observe() const {
        return _kernel->observe();
    }

    /**
     * Registers a predefined observable with the kernel, which enables evaluation during the simulation.
     * @param observable the observable instance
     * @return an observable handle that allows for post-hoc modification of the observable
     */
    ObservableHandle registerObservable(std::unique_ptr<readdy::model::observables::ObservableBase> observable) {
        return _kernel->registerObservable(std::move(observable));
    }

    /**
     * A method to access the particle positions of a certain type.
     * @param type the type
     * @return a vector containing the particle positions
     */
    [[nodiscard]] std::vector<Vec3> getParticlePositions(const std::string &type) {
        auto typeId = _kernel->context().particleTypes().idOf(type);
        const auto particles = _kernel->stateModel().getParticles();
        std::vector<Vec3> positions;
        positions.reserve(particles.size());
        for (const auto &p : particles) {
            if (p.type() == typeId) {
                positions.push_back(p.pos());
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
            readdy::model::Particle p{x, y, z, context().particleTypes().idOf(type)};
            _kernel->actions().addParticles(p)->perform();
        } else {
            log::error("particle position was not in bounds of the simulation box!");
        }
    }

    /**
     * Method that gives access to all the positions of all the particles in the system.
     * @return a vector containing all the positions
     */
    [[nodiscard]] std::vector<Vec3> getAllParticlePositions() const {
        return _kernel->stateModel().getParticlePositions();
    }

    /**
     * Yields a nonmodifiable reference to the current context of this simulation.
     * @return cref to the context
     */
    [[nodiscard]] const model::Context &context() const {
        return _kernel->context();
    }

    /**
     * Method that returns the type of the currently selected kernel.
     * Raises if there is no kernel selected.
     * @return the kernel id string
     */
    [[nodiscard]] const std::string &selectedKernelType() const {
        return _kernel->name();
    }

    /**
     * Runs the simulation using the default simulation loop for a given number of time steps and a given time step.
     * @param steps the number of steps
     * @param timeStep the time step
     * @see runScheme()
     */
    virtual void run(TimeStep steps, scalar timeStep) {
        createLoop(timeStep).run(steps);
    }

    api::SimulationLoop createLoop(scalar timeStep) {
        return api::SimulationLoop(_kernel.get(), timeStep);
    }

    /**
     * Checks if the kernel is running on single precision.
     * @return true if it is running on single precision
     */
    [[nodiscard]] bool singlePrecision() const {
        return _kernel->singlePrecision();
    }

    /**
     * Checks if the kernel is running on double precision.
     * @return true if it is running on double precision
     */
    [[nodiscard]] bool doublePrecision() const {
        return _kernel->doublePrecision();
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

    [[nodiscard]] const model::StateModel &stateModel() const {
        return _kernel->stateModel();
    }

    readdy::model::actions::ActionFactory &actions() {
        return _kernel->actions(); // RETURRRRRN
    }

private:
    plugin::KernelProvider::kernel_ptr _kernel;
};

}
