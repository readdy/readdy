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
                        using _reaction_idx_t = _event_t::index_type;

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

                    class GillespieParallel : public readdy::model::programs::reactions::GillespieParallel {
                        using kernel_t = readdy::kernel::cpu::CPUKernel;
                        using vec_t = readdy::model::Vec3;
                        using data_t = decltype(std::declval<kernel_t>().getKernelStateModel().getParticleData());
                        using nl_t = decltype(std::declval<kernel_t>().getKernelStateModel().getNeighborList());
                        using ctx_t = std::remove_const<decltype(std::declval<kernel_t>().getKernelContext())>::type;
                    public:
                        GillespieParallel(kernel_t const *const kernel);
                        ~GillespieParallel();
                        virtual void execute() override;
                        void clear();
                    private:
                        kernel_t const *const kernel;
                        double maxReactionRadius = 0.0;
                        double boxWidth = 0.0;
                        unsigned int longestAxis;
                        unsigned int otherAxis1, otherAxis2;
                        struct HaloBox;
                        std::vector<HaloBox> boxes;
                        std::vector<unsigned long> problematicParticles {};

                        /**
                         * look for the longest axis and divide it into n_threads parts, yielding halo boxes.
                         */
                        void setupBoxes();
                        /**
                         * Sorts particles into boxes and adds them to the problematic particles should they fall into
                         * a halo region
                         */
                        void fillBoxes();
                        /**
                         * Executes the gillespie algorithm for each box and gives an update on problematic particles
                         */
                        void handleBoxReactions();

                        void handleProblematic(const unsigned long idx, const HaloBox &box, const ctx_t ctx, data_t data, nl_t nl, std::vector<unsigned long>& update) const;

                    };
                }

            }
        }
    }
}
#endif //READDY_MAIN_REACTIONS_H_H
