/**
 * << detailed description >>
 *
 * @file CPUTestParallelizedGillespie.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 01.09.16
 */

#include <gtest/gtest.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace {

    TEST(TestParallelGillespie, Sanity) {
        readdy::kernel::cpu::CPUKernel kernel;
        kernel.getKernelContext().setBoxSize(10, 10, 11);
        kernel.getKernelContext().setDiffusionConstant("A", 10.0);
        kernel.registerReaction<readdy::model::reactions::Fusion>("Fusion", "A", "A", "A", 10, 1.0);
        kernel.addParticle("A", {-5, .2, -5.5});
        kernel.addParticle("A", {-5, .2, 5.5});
        kernel.addParticle("A", {-5, .2, 0});
        kernel.getKernelContext().configure();
        auto prog = kernel.createProgram<readdy::model::programs::reactions::GillespieParallel>();
        prog->execute();
    }
}