
/**
 * << detailed description >>
 *
 * @file CPUMiscActions.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 17.03.20
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/api/Saver.h>
#include "../CPUKernel.h"

namespace readdy::kernel::cpu::actions {
class CPUEvaluateObservables : public readdy::model::actions::EvaluateObservables {
public:
    explicit CPUEvaluateObservables(CPUKernel *kernel) : kernel(kernel) {}

    void perform(TimeStep t) override {
        kernel->evaluateObservables(t);
    }

private:
    CPUKernel *kernel;
};

class CPUMakeCheckpoint : public readdy::model::actions::MakeCheckpoint {
public:
    CPUMakeCheckpoint(CPUKernel *kernel, const std::string& base, std::size_t maxNSaves) : kernel(kernel), saver(base, maxNSaves) {}

    void perform(TimeStep t) override {
        saver.makeCheckpoint(kernel, t);
    }
private:
    CPUKernel *kernel;
    readdy::api::Saver saver;
};

}
