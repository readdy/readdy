/**
 * << detailed description >>
 *
 * @file SimulationScheme.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.08.16
 */

#ifndef READDY_MAIN_SIMULATIONSCHEME_H
#define READDY_MAIN_SIMULATIONSCHEME_H

#include <memory>
#include <readdy/common/Types.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/programs/Programs.h>

namespace readdy {
    namespace api {
        namespace {
            struct SimulationScheme {
                SimulationScheme(model::Kernel *const kernel) : kernel(kernel) {}

                virtual void run(const model::time_step_type steps) = 0;

            protected:
                model::Kernel *const kernel;
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

            private:
                friend class SchemeBuilder<ReaDDyScheme>;
                std::unique_ptr<model::programs::Program> integrator = nullptr;
                std::unique_ptr<model::programs::Program> forces = nullptr;
                std::unique_ptr<model::programs::Program> reactionScheduler = nullptr;
                std::unique_ptr<model::programs::UpdateNeighborList> neighborList = nullptr;
                bool evaluateObservables = true;
            };

            template<typename SchemeType>
            class SchemeBuilder {
            public:

                SchemeBuilder(model::Kernel *const kernel, bool useDefaults = true) : scheme(SchemeType(kernel)),
                                                                                      useDefaults(useDefaults) {}

                SchemeBuilder &withIntegrator(std::unique_ptr<model::programs::Program> integrator) {
                    scheme.integrator = std::move(integrator);
                    return *this;
                }

                SchemeBuilder &withReactionScheduler(std::unique_ptr<model::programs::Program> reactionScheduler) {
                    scheme.reactionScheduler = std::move(reactionScheduler);
                    return *this;
                }

                SchemeBuilder &evaluateObservables(bool evaluate) {
                    scheme.evaluateObservables = evaluate;
                    evaluateObservablesSet = true;
                    return *this;
                }

                SchemeBuilder &includeForces() {
                    scheme.forces = scheme.kernel->createProgram<readdy::model::programs::CalculateForces>();
                    return *this;
                }

                SchemeType build() {
                    if (useDefaults) {
                        if (!scheme.integrator) {
                            withIntegrator(scheme.kernel->createProgram<readdy::model::programs::EulerBDIntegrator>());
                        }
                        if (!scheme.reactionScheduler) {
                            withReactionScheduler(
                                    scheme.kernel->createProgram<readdy::model::programs::reactions::Gillespie>());
                        }
                        if (!evaluateObservablesSet) {
                            evaluateObservables(true);
                        }
                        if (!scheme.forces) {
                            includeForces();
                        }
                    }
                    if (scheme.forces || scheme.reactionScheduler) {
                        scheme.neighborList = scheme.kernel->createProgram<readdy::model::programs::UpdateNeighborList>();
                    }
                    return std::move(scheme);
                }

            private:
                const bool useDefaults;
                bool evaluateObservablesSet = false;
                SchemeType scheme;
            };
        }
    }
}
#endif //READDY_MAIN_SIMULATIONSCHEME_H
