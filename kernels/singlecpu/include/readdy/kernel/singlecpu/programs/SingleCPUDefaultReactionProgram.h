/**
 * << detailed description >>
 *
 * @file SingleCPUDefaultReactionProgram.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#ifndef READDY_MAIN_SINGLECPUDEFAULTREACTIONPROGRAM_H
#define READDY_MAIN_SINGLECPUDEFAULTREACTIONPROGRAM_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                class SingleCPUDefaultReactionProgram : public readdy::model::programs::DefaultReactionProgram {

                public:
                    SingleCPUDefaultReactionProgram(SingleCPUKernel const *const kernel);

                    virtual void execute() override;

                    virtual void registerReactionScheme_11(const std::string &reactionName, reaction_11 fun) override;
                    virtual void registerReactionScheme_12(const std::string &reactionName, reaction_12 fun) override;
                    virtual void registerReactionScheme_21(const std::string &reactionName, reaction_21 fun) override;
                    virtual void registerReactionScheme_22(const std::string &reactionName, reaction_22 fun) override;

                protected:
                    SingleCPUKernel const *const kernel;
                    std::unordered_map<boost::uuids::uuid, reaction_11, boost::hash<boost::uuids::uuid>> mapping_11 {};
                    std::unordered_map<boost::uuids::uuid, reaction_12, boost::hash<boost::uuids::uuid>> mapping_12 {};
                    std::unordered_map<boost::uuids::uuid, reaction_21, boost::hash<boost::uuids::uuid>> mapping_21 {};
                    std::unordered_map<boost::uuids::uuid, reaction_22, boost::hash<boost::uuids::uuid>> mapping_22 {};
                };
            }
        }
    }
}
#endif //READDY_MAIN_SINGLECPUDEFAULTREACTIONPROGRAM_H
