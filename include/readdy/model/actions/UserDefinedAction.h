//
// Created by mho on 7/2/18.
//

#ifndef READDY_MAIN_USERDEFINEDACTION_H
#define READDY_MAIN_USERDEFINEDACTION_H

#include <readdy/model/Kernel.h>
#include "Action.h"

namespace readdy {
namespace model {
namespace actions {

class UserDefinedAction : public TimeStepDependentAction {
public:
    explicit UserDefinedAction(scalar timeStep) : TimeStepDependentAction(timeStep) {}

    model::Kernel*& kernel() {
        return _kernel;
    }
    
    model::Kernel* const& kernel() const {
        return _kernel;
    }
    
protected:
    model::Kernel *_kernel {nullptr};
};

}
}
}

#endif //READDY_MAIN_USERDEFINEDACTION_H
