/**
 * Verify that the programs of the different kernels give the same results,
 * e.g. that calculated forces and possible reaction events are identical. These
 * tests work on the level of programs and thus do not check for kernel-specific
 * implementation details.
 *
 * @file TestCompareKernels.cpp
 * @brief Assure semantic correctness and equality of kernel implementations/programs.
 * @author chrisfroe
 * @date 08.08.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>

namespace {
    /**
     * Small system, partly periodic, different x,y,z boxlength.
     * particle radius 0.9 and harmonic repulsion (kappa=1) -> maxcutoff = 1.8
     * simbox = 2x4x6 -> numberBoxes = 1x2x3, periodic = false x true x true
     * particle positions:
     * [(0.01,1.99,0.01), (1.99,2.01,0.01), (0.5,0.01,2.01), (0.5,3.99,2.01), (0.01,2.01,4.01), (1.99,1.99,4.01), (0.01,1.99,5.99)]
     * force/energy contributions from particles (2,3) and (6,0) both with distance 0.02
     * f2=(0,1.78,0) f3=(0,-1.78,0) f6=(0,0,-1.78) f0=(0,0,1.78)
     * e=0.5 * 1.78**2 , totalenergy=2*e
     */
    TEST(TestCompareKernels, ForcesSmallSystem) {
        auto singleCpuKernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        auto cpuKernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    }

    TEST(TestCompareKernels, ForcesLargerSystem) {
        auto singleCpuKernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        auto cpuKernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    }
}
