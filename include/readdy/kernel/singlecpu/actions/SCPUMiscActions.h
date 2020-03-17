
/**
 * << detailed description >>
 *
 * @file SCPUMiscActions.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 17.03.20
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/api/Saver.h>

namespace readdy::kernel::scpu::actions {
class SCPUEvaluateObservables : public readdy::model::actions::EvaluateObservables {
public:
    explicit SCPUEvaluateObservables(SCPUKernel *kernel) : kernel(kernel) {}

    void perform(TimeStep t) override {
        kernel->evaluateObservables(t);
    }

private:
    SCPUKernel *kernel;
};

class SCPUMakeCheckpoint : public readdy::model::actions::MakeCheckpoint {
public:
    SCPUMakeCheckpoint(SCPUKernel *kernel, const std::string& base, std::size_t maxNSaves) : kernel(kernel), saver(base, maxNSaves) {}

    void perform(TimeStep t) override {
        saver.makeCheckpoint(kernel, t);
    }
private:
    SCPUKernel *kernel;
    readdy::api::Saver saver;
};

}