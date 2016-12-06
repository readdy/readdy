//
// Created by clonker on 11.04.16.
//

#include <readdy/kernel/singlecpu/programs/SCPUTestProgram.h>
#include <readdy/common/logging.h>

namespace sctp = readdy::kernel::scpu::programs;

sctp::SCPUTestProgram::SCPUTestProgram() : readdy::model::programs::Test() {

}

void sctp::SCPUTestProgram::execute() {
    log::console()->debug("execute called!");
}

sctp::SCPUTestProgram::~SCPUTestProgram() = default;







