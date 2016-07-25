/**
 * << detailed description >>
 *
 * @file KernelTest.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.07.16
 */


#ifndef READDY_MAIN_KERNELTEST_H
#define READDY_MAIN_KERNELTEST_H

#include <gtest/gtest.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>

class KernelTest : public ::testing::TestWithParam<std::string> {
public:
    const std::unique_ptr<readdy::model::Kernel> kernel;

    explicit KernelTest() : kernel(readdy::plugin::KernelProvider::getInstance().create(GetParam())){
    }
};

#endif //READDY_MAIN_KERNELTEST_H
