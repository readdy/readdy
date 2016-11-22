/**
 * << detailed description >>
 *
 * @file UncontrolledApproximation.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_DENSE_UNCONTROLLEDAPPROXIMATION_H
#define READDY_DENSE_UNCONTROLLEDAPPROXIMATION_H

#include <readdy/kernel/cpu_dense/Kernel.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
namespace reactions {

class UncontrolledApproximation : public readdy::model::programs::reactions::UncontrolledApproximation {
public:
    UncontrolledApproximation(const Kernel *const kernel);

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
    Kernel const *const kernel;
};
}
}
}
}
}
#endif //READDY_DENSE_UNCONTROLLEDAPPROXIMATION_H
