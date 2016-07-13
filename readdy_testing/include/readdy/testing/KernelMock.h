/**
 * << detailed description >>
 *
 * @file KernelMock.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_MAIN_KERNELMOCK_H
#define READDY_MAIN_KERNELMOCK_H

#include <readdy/model/Kernel.h>
#include <gmock/gmock.h>

namespace readdy {
    namespace testing {
        class KernelMock : public readdy::model::Kernel {

        public:
            KernelMock(const std::string &name) : Kernel(name) { }

            MOCK_CONST_METHOD0(getProgramFactory, readdy::model::programs::ProgramFactory&(void));
            MOCK_CONST_METHOD0(getKernelStateModel, readdy::model::KernelStateModel&(void));
            MOCK_CONST_METHOD0(getKernelContext, readdy::model::KernelContext&(void));
            MOCK_CONST_METHOD0(getPotentialFactory, readdy::model::potentials::PotentialFactory&(void));
            MOCK_CONST_METHOD0(getReactionFactory, readdy::model::reactions::ReactionFactory&(void));

        };
    }
}
#endif //READDY_MAIN_KERNELMOCK_H
