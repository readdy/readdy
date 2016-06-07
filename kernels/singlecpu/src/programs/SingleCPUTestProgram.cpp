//
// Created by clonker on 11.04.16.
//

#include <boost/log/trivial.hpp>
#include <readdy/kernel/singlecpu/programs/SingleCPUTestProgram.h>

namespace sctp = readdy::kernel::singlecpu::programs;

sctp::SingleCPUTestProgram::SingleCPUTestProgram() : readdy::model::TestProgram() {

}

void sctp::SingleCPUTestProgram::execute() {
    BOOST_LOG_TRIVIAL(debug) << "execute called!";
}

sctp::SingleCPUTestProgram::~SingleCPUTestProgram() = default;







