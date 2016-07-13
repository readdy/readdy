/**
 * << detailed description >>
 *
 * @file CPUDefaultReactionProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#ifndef READDY_MAIN_CPUDEFAULTREACTIONPROGRAM_H
#define READDY_MAIN_CPUDEFAULTREACTIONPROGRAM_H

#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/model/programs/Programs.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                namespace reactions {
                    class UncontrolledApproximation : public readdy::model::programs::reactions::UncontrolledApproximation {
                    public:
                        UncontrolledApproximation(const CPUKernel *const kernel);

                        virtual void execute() override;

                        virtual void registerReactionScheme_11(const std::string &reactionName,
                                                               reaction_11 fun) override {

                        }

                        virtual void registerReactionScheme_12(const std::string &reactionName,
                                                               reaction_12 fun) override {

                        }

                        virtual void registerReactionScheme_21(const std::string &reactionName,
                                                               reaction_21 fun) override {

                        }

                        virtual void registerReactionScheme_22(const std::string &reactionName,
                                                               reaction_22 fun) override {

                        }


                    protected:
                        CPUKernel const *const kernel;
                    };
                }

            }
        }
    }
}
#endif //READDY_MAIN_CPUDEFAULTREACTIONPROGRAM_H
