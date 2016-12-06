/**
 * << detailed description >>
 *
 * @file CPUProgramFactory.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */


#ifndef READDY_DENSE_CPUPROGRAMFACTORY_H
#define READDY_DENSE_CPUPROGRAMFACTORY_H

#include <readdy/kernel/cpu_dense/CPUDKernel.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
class CPUDProgramFactory : public readdy::model::programs::ProgramFactory {

public:
    CPUDProgramFactory(CPUDKernel *kernel);
};

}
}
}
}

#endif //READDY_DENSE_CPUPROGRAMFACTORY_H
