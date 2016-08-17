/**
 * << detailed description >>
 *
 * @file KernelContextMock.h
 * @brief Mock a KernelContext if only a few functions of it are required without creating a whole kernel.
 * @author chrisfroe
 * @date 17.08.16
 */

#ifndef READDY_MAIN_KERNELCONTEXTMOCK_H
#define READDY_MAIN_KERNELCONTEXTMOCK_H

#include <gmock/gmock.h>

namespace readdy {
    namespace testing {
        /**
         * Since the KernelContext has no virtual methods, it cannot be simply mocked. Instead create another
         * class with the test-relevant functions only. Then the function/class to be tested must be templated with
         * respect to the KernelContext type. Note that this class is not related to the class to-be-mocked in any way.
         * This implies that the dev'er is responsible for keeping these methods congruent with the "real" code.
         * See https://github.com/google/googletest/blob/master/googlemock/docs/CookBook.md#mocking-nonvirtual-methods
         */
        class KernelContextMock {

        public:
            MOCK_CONST_METHOD0(getBoxSize, std::array<double,3>&(void));
            MOCK_CONST_METHOD0(getVectorAllOrder2Potentials, const std::vector<readdy::model::potentials::PotentialOrder2*>(void));
            MOCK_CONST_METHOD0(getAllOrder2Reactions, const std::vector<const readdy::model::reactions::Reaction<2> *>(void));
            MOCK_CONST_METHOD0(getPeriodicBoundary, const std::array<bool, 3>&(void));
        };
    }
}

#endif //READDY_MAIN_KERNELCONTEXTMOCK_H
