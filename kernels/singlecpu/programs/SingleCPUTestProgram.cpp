//
// Created by clonker on 11.04.16.
//

#include <boost/log/trivial.hpp>
#include "SingleCPUTestProgram.h"

using sctp =readdy::kernel::singlecpu::programs::SingleCPUTestProgram;

sctp::SingleCPUTestProgram() : readdy::plugin::TestProgram(){

}

void sctp::execute() {
    BOOST_LOG_TRIVIAL(debug) << "execute called!";
}





