/**
 * << detailed description >>
 *
 * @file Kernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#include "readdy/kernel/cpu_dense/Kernel.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {

const std::string Kernel::name = "CPU Dense";

}
}
}


const char* name()  {
    return readdy::kernel::cpu_dense::Kernel::name.c_str();
}

readdy::model::Kernel* createKernel() {
    return readdy::kernel::cpu_dense::Kernel::create();
}
