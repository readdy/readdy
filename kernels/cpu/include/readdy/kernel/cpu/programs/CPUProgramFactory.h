/**
 * << detailed description >>
 *
 * @file CPUProgramFactory.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */


#ifndef READDY_MAIN_CPUPROGRAMFACTORY_H
#define READDY_MAIN_CPUPROGRAMFACTORY_H

#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                class CPUProgramFactory : public readdy::model::programs::ProgramFactory {

                public:
                    CPUProgramFactory(CPUKernel *kernel);
                };

            }
        }
    }
}

#endif //READDY_MAIN_CPUPROGRAMFACTORY_H
