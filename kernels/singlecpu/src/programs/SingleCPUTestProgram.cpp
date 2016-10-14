//
// Created by clonker on 11.04.16.
//

#include <readdy/kernel/singlecpu/programs/SingleCPUTestProgram.h>
#include <readdy/common/logging.h>

namespace sctp = readdy::kernel::singlecpu::programs;

sctp::SingleCPUTestProgram::SingleCPUTestProgram() : readdy::model::programs::Test() {

}

void sctp::SingleCPUTestProgram::execute() {
    log::console()->debug("execute called!");
}

sctp::SingleCPUTestProgram::~SingleCPUTestProgram() = default;







