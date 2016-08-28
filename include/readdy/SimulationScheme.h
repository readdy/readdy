/**
 * Vital parts of the scheme api are defined in this header. First of all, there is the SimulationScheme. Its task is
 * to execute an integrator(program), a force(program), a reaction scheduler(program) and a neighborList (program)
 * in a certain way. One example would be the ReaDDyScheme, which will configure the context, then:
 *
 * - create neighbor list
 * - calcualte forces
 * - evaluate observables (0)
 * - for i in range(steps):
 *      - evaluate integrator
 *      - update neighbor list
 *      - calculate forces
 *      - evaluate reactions
 *      - update the neighbor list
 *      - calcualte forces
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
#include <readdy/common/Types.h>
#include <readdy/model/Kernel.h>

namespace readdy {
    namespace api {
        template<typename SchemeType>
        class SchemeConfigurator;

        struct SimulationScheme {
            SimulationScheme(model::Kernel *const kernel) : kernel(kernel) {}

            virtual void run(const model::time_step_type steps) = 0;

        protected:
            template<typename SchemeType>
            friend class SchemeConfigurator;

            model::Kernel *const kernel;
            std::unique_ptr<model::programs::Program> integrator = nullptr;
            std::unique_ptr<model::programs::Program> forces = nullptr;
            std::unique_ptr<model::programs::Program> reactionScheduler = nullptr;
            std::unique_ptr<model::programs::UpdateNeighborList> neighborList = nullptr;
            bool evaluateObservables = true;
        };

        class ReaDDyScheme : public SimulationScheme {
        public:
            ReaDDyScheme(model::Kernel *const kernel) : SimulationScheme(kernel) {};

            virtual void run(const model::time_step_type steps) override {
                kernel->getKernelContext().configure();

                if (neighborList) neighborList->execute();
                if (forces) forces->execute();
                if (evaluateObservables) kernel->evaluateObservables(0);
                for (model::time_step_type &&t = 0; t < steps; ++t) {
                    if (integrator) integrator->execute();
                    if (neighborList) neighborList->execute();
                    if (forces) forces->execute();

                    if (reactionScheduler) reactionScheduler->execute();
                    if (neighborList) neighborList->execute();
                    if (forces) forces->execute();
                    if (evaluateObservables) kernel->evaluateObservables(t + 1);
                }

                if (neighborList) neighborList->setAction(model::programs::UpdateNeighborList::Action::clear);
                if (neighborList) neighborList->execute();
            }
        };

        template<typename SchemeType>
        class SchemeConfigurator {
            static_assert(std::is_base_of<SimulationScheme, SchemeType>::value,
                          "SchemeType must inherit from readdy::api::SimulationScheme");
        public:

            SchemeConfigurator(model::Kernel *const kernel, bool useDefaults = true) : scheme(std::make_unique<SchemeType>(kernel)),
                                                                                       useDefaults(useDefaults) {}
            SchemeConfigurator& withIntegrator(std::unique_ptr<model::programs::Program> integrator) {
                scheme->integrator = std::move(integrator);
                return *this;
            }

            template<typename IntegratorType>
            SchemeConfigurator& withIntegrator() {
                scheme->integrator = scheme->kernel->template createProgram<IntegratorType>();
                return *this;
            }

            SchemeConfigurator& withIntegrator(const std::string &integratorName) {
                scheme->integrator = scheme->kernel->createProgram(integratorName);
                return *this;
            }

            template<typename ReactionSchedulerType>
            SchemeConfigurator& withReactionScheduler() {
                scheme->reactionScheduler = scheme->kernel->template createProgram<ReactionSchedulerType>();
                return *this;
            }

            SchemeConfigurator& withReactionScheduler(std::unique_ptr<model::programs::Program> reactionScheduler) {
                scheme->reactionScheduler = std::move(reactionScheduler);
                return *this;
            }

            SchemeConfigurator& withReactionScheduler(const std::string &schedulerName) {
                scheme->reactionScheduler = scheme->kernel->createProgram(schedulerName);
                return *this;
            }

            SchemeConfigurator& evaluateObservables(bool evaluate = true) {
                scheme->evaluateObservables = evaluate;
                evaluateObservablesSet = true;
                return *this;
            }

            SchemeConfigurator& includeForces(bool include = true) {
                if (include) {
                    scheme->forces = scheme->kernel->template createProgram<readdy::model::programs::CalculateForces>();
                } else {
                    scheme->forces = nullptr;
                }
                includeForcesSet = true;
                return *this;
            }

            std::unique_ptr<SchemeType> configure() {
                if (useDefaults) {
                    if (scheme->integrator) {
                        scheme->integrator = scheme->kernel->template createProgram<readdy::model::programs::EulerBDIntegrator>();
                    }
                    if (!scheme->reactionScheduler) {
                        scheme->reactionScheduler = scheme->kernel->template createProgram<readdy::model::programs::reactions::Gillespie>();
                    }
                    if (!evaluateObservablesSet) {
                        scheme->evaluateObservables = true;
                    }
                    if (!scheme->forces && !includeForcesSet) {
                        scheme->forces = scheme->kernel->template createProgram<readdy::model::programs::CalculateForces>();
                    }
                }
                if (scheme->forces || scheme->reactionScheduler) {
                    scheme->neighborList = scheme->kernel
                            ->template createProgram<readdy::model::programs::UpdateNeighborList>();
                }
                std::unique_ptr<SchemeType> ptr = std::move(scheme);
                scheme = nullptr;
                return ptr;
            }

            void configureAndRun(const model::time_step_type steps) {
                configure()->run(steps);
            }

        private:
            const bool useDefaults;
            bool evaluateObservablesSet = false;
            bool includeForcesSet = false;
            std::unique_ptr<SchemeType> scheme = nullptr;
        };
    }
}
#endif //READDY_MAIN_SIMULATIONSCHEME_H
