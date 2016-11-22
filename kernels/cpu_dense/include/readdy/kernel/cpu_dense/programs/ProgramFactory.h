/**
 * << detailed description >>
 *
 * @file CPUProgramFactory.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */


#ifndef READDY_DENSE_CPUPROGRAMFACTORY_H
#define READDY_DENSE_CPUPROGRAMFACTORY_H

#include <readdy/kernel/cpu_dense/Kernel.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
class ProgramFactory : public readdy::model::programs::ProgramFactory {

public:
    ProgramFactory(Kernel *kernel);
};

}
}
}
}

#endif //READDY_DENSE_CPUPROGRAMFACTORY_H
