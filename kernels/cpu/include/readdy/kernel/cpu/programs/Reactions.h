/**
 * << detailed description >>
 *
 * @file Reactions.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.07.16
 */

#ifndef READDY_MAIN_REACTIONS_H_H
#define READDY_MAIN_REACTIONS_H_H
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                namespace reactions {
                    class UncontrolledApproximation : public readdy::model::programs::reactions::UncontrolledApproximation {
                    public:
                        UncontrolledApproximation(const CPUKernel *const kernel);

                        virtual void execute() override;
                        virtual void registerReactionScheme_11(const std::string &reactionName, reaction_11 fun) override {
                            throw std::runtime_error("not supported for cpu kernel thus far");
                        }

                        virtual void registerReactionScheme_12(const std::string &reactionName, reaction_12 fun) override {
                            throw std::runtime_error("not supported for cpu kernel thus far");
                        }
                        virtual void registerReactionScheme_21(const std::string &reactionName, reaction_21 fun) override {
                            throw std::runtime_error("not supported for cpu kernel thus far");

                        }
                        virtual void registerReactionScheme_22(const std::string &reactionName, reaction_22 fun) override {
                            throw std::runtime_error("not supported for cpu kernel thus far");
                        }
                    protected:
                        CPUKernel const *const kernel;
                    };

                    class Gillespie : public readdy::model::programs::reactions::Gillespie {
                        using _event_t = readdy::kernel::singlecpu::programs::reactions::ReactionEvent;
                        using _reaction_idx_t = _event_t ::index_type;
                    public:

                        Gillespie(CPUKernel const *const kernel);

                        virtual void execute() override {
                            const auto& ctx = kernel->getKernelContext();
                            auto data = kernel->getKernelStateModel().getParticleData();
                            const auto &dist = ctx.getDistSquaredFun();
                            const auto &fixPos = ctx.getFixPositionFun();

                            double alpha = 0.0;
                            auto events = gatherEvents(alpha);
                            auto newParticles = handleEvents(std::move(events), alpha);

                            // reposition particles to respect the periodic b.c.
                            std::for_each(newParticles.begin(), newParticles.end(),
                                          [&fixPos](readdy::model::Particle &p) { fixPos(p.getPos()); });

                            // update data structure
                            data->deactivateMarked();
                            data->addParticles(newParticles);
                        }

                    protected:
                        virtual std::vector<_event_t> gatherEvents(double& alpha);
                        virtual std::vector<readdy::model::Particle> handleEvents(std::vector<_event_t> events, double alpha);
                        CPUKernel const *const kernel;
                    };
                }

            }
        }
    }
}
#endif //READDY_MAIN_REACTIONS_H_H
