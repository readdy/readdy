/**
 * << detailed description >>
 *
 * @file UncontrolledApproximation.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_CPUKERNEL_UNCONTROLLEDAPPROXIMATION_H
#define READDY_CPUKERNEL_UNCONTROLLEDAPPROXIMATION_H

#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

class CPUUncontrolledApproximation : public readdy::model::programs::reactions::UncontrolledApproximation {
public:
    CPUUncontrolledApproximation(const CPUKernel *const kernel);

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
}
}
}
}
}
#endif //READDY_CPUKERNEL_UNCONTROLLEDAPPROXIMATION_H
