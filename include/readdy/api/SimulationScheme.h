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
 * Vital parts of the scheme api are defined in this header. First of all, there is the SimulationScheme. Its task is
 * to execute an integrator(program), a force(program), a reaction scheduler(program) and a neighborList (program)
 * in a certain way. One example would be the ReaDDyScheme, which will configure the context, then:
 *
 * - create neighbor list
 * - calculate forces
 * - evaluate observables (0)
 * - for i in range(steps):
 *      - evaluate integrator
 *      - update neighbor list
 *      - calculate forces
 *      - evaluate reactions
 *      - update the neighbor list
 *      - calculate forces
 *      - evaluate observables (i+1)
 *  - clear neighbor list
 *
 *  The configurator makes use of the builder pattern and either sets default values or not, depending on
 *  its constructor argument.
 *
 * @file SimulationScheme.h
 * @brief Header file containing the SchemeConfigurator<T> and a corresponding SimulationScheme superclass definition.
 * @author clonker
 * @date 23.08.16
 */

#ifndef READDY_MAIN_SIMULATIONSCHEME_H
#define READDY_MAIN_SIMULATIONSCHEME_H

#include <memory>
#include <type_traits>
#include <readdy/model/Kernel.h>

namespace readdy {
    namespace api {
        template<typename SchemeType>
        class SchemeConfigurator;

        struct SimulationScheme {
            SimulationScheme(model::Kernel *const kernel) : kernel(kernel) {}

            virtual void run(const model::observables::time_step_type steps) = 0;

        protected:
            template<typename SchemeType>
            friend class SchemeConfigurator;

            model::Kernel *const kernel;
            std::unique_ptr<model::actions::TimeStepDependentAction> integrator = nullptr;
            std::unique_ptr<model::actions::Action> forces = nullptr;
            std::unique_ptr<model::actions::TimeStepDependentAction> reactionScheduler = nullptr;
            std::unique_ptr<model::actions::UpdateNeighborList> neighborList = nullptr;
            std::unique_ptr<model::actions::UpdateNeighborList> clearNeighborList = nullptr;
            bool evaluateObservables = true;
            model::observables::time_step_type start = 0;
        };

        class ReaDDyScheme : public SimulationScheme {
        public:
            ReaDDyScheme(model::Kernel *const kernel) : SimulationScheme(kernel) {};

            virtual void run(const model::observables::time_step_type steps) override {
                kernel->getKernelContext().configure(true);

                if (neighborList) neighborList->perform();
                if (forces) forces->perform();
                if (evaluateObservables) kernel->evaluateObservables(start);
                // show every 1% of the simulation
                const std::size_t progressOutputStride = static_cast<std::size_t>(steps / 100);
                for (model::observables::time_step_type t = start; t < start+steps; ++t) {
                    if (integrator) integrator->perform();
                    if (neighborList) neighborList->perform();
                    // if (forces) forces->perform();

                    if (reactionScheduler) reactionScheduler->perform();
                    if (neighborList) neighborList->perform();
                    if (forces) forces->perform();
                    if (evaluateObservables) kernel->evaluateObservables(t + 1);
                    if(progressOutputStride > 0 && (t - start) % progressOutputStride == 0) {
                        log::console()->debug("Simulation progress: {} / {} steps", (t - start), steps);
                    }
                }

                if (clearNeighborList) clearNeighborList->perform();
                start += steps;
                log::console()->debug("Simulation completed");
            }
        };

        template<typename SchemeType>
        class SchemeConfigurator {
            static_assert(std::is_base_of<SimulationScheme, SchemeType>::value,
                          "SchemeType must inherit from readdy::api::SimulationScheme");
        public:

            SchemeConfigurator(model::Kernel *const kernel, bool useDefaults = true) : scheme(std::make_unique<SchemeType>(kernel)),
                                                                                       useDefaults(useDefaults) {}
            SchemeConfigurator& withIntegrator(std::unique_ptr<model::actions::TimeStepDependentAction> integrator) {
                scheme->integrator = std::move(integrator);
                return *this;
            }

            template<typename IntegratorType>
            SchemeConfigurator& withIntegrator() {
                scheme->integrator = scheme->kernel->template createAction<IntegratorType>(0.);
                return *this;
            }

            SchemeConfigurator& withIntegrator(const std::string &integratorName) {
                scheme->integrator = scheme->kernel->getActionFactory().createIntegrator(integratorName, 0.);
                return *this;
            }

            template<typename ReactionSchedulerType>
            SchemeConfigurator& withReactionScheduler() {
                scheme->reactionScheduler = scheme->kernel->template createAction<ReactionSchedulerType>(0.);
                return *this;
            }

            SchemeConfigurator& withReactionScheduler(std::unique_ptr<model::actions::TimeStepDependentAction> reactionScheduler) {
                scheme->reactionScheduler = std::move(reactionScheduler);
                return *this;
            }

            SchemeConfigurator& withReactionScheduler(const std::string &name) {
                scheme->reactionScheduler = scheme->kernel->getActionFactory().createReactionScheduler(name, 0.);
                return *this;
            }

            SchemeConfigurator& evaluateObservables(bool evaluate = true) {
                scheme->evaluateObservables = evaluate;
                evaluateObservablesSet = true;
                return *this;
            }

            SchemeConfigurator& includeForces(bool include = true) {
                if (include) {
                    scheme->forces = scheme->kernel->template createAction<readdy::model::actions::CalculateForces>();
                } else {
                    scheme->forces = nullptr;
                }
                includeForcesSet = true;
                return *this;
            }

            std::unique_ptr<SchemeType> configure(double timeStep) {
                using default_integrator_t = readdy::model::actions::EulerBDIntegrator;
                using default_reactions_t = readdy::model::actions::reactions::Gillespie;
                using calculate_forces_t = readdy::model::actions::CalculateForces;
                using update_neighbor_list_t = readdy::model::actions::UpdateNeighborList;
                if (useDefaults) {
                    if (!scheme->integrator) {
                        scheme->integrator = scheme->kernel->template createAction<default_integrator_t>(timeStep);
                    }
                    if (!scheme->reactionScheduler) {
                        scheme->reactionScheduler = scheme->kernel->template createAction<default_reactions_t>(timeStep);
                    }
                    if (!evaluateObservablesSet) {
                        scheme->evaluateObservables = true;
                    }
                    if (!scheme->forces && !includeForcesSet) {
                        scheme->forces = scheme->kernel->template createAction<calculate_forces_t>();
                    }
                }
                if (scheme->forces || scheme->reactionScheduler) {
                    scheme->neighborList = scheme->kernel
                            ->template createAction<update_neighbor_list_t>(update_neighbor_list_t::Operation::create, -1);
                    scheme->clearNeighborList = scheme->kernel
                            ->template createAction<update_neighbor_list_t>(update_neighbor_list_t::Operation::clear, -1);
                }
                if(scheme->integrator) scheme->integrator->setTimeStep(timeStep);
                if(scheme->reactionScheduler) scheme->reactionScheduler->setTimeStep(timeStep);
                std::unique_ptr<SchemeType> ptr = std::move(scheme);
                scheme = nullptr;
                return ptr;
            }

            void configureAndRun(double timeStep, const model::observables::time_step_type steps) {
                configure(timeStep)->run(steps);
            }

        private:
            bool useDefaults;
            bool evaluateObservablesSet = false;
            bool includeForcesSet = false;
            std::unique_ptr<SchemeType> scheme = nullptr;
        };
    }
}
#endif //READDY_MAIN_SIMULATIONSCHEME_H
