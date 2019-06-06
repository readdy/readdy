/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * @file TestBreakingBonds.cpp
 * @brief Test implementation-independent execution of action that breaks bonds
 * @author chrisfroe
 * @date 27.05.19
 */

#include <catch2/catch.hpp>

#include <readdy/testing/KernelTest.h>

namespace m = readdy::model;
using namespace readdytesting::kernel;

TEMPLATE_TEST_CASE("Test breaking bonds.", "[breakbonds]", SingleCPU, CPU) {
    auto kernel = readdytesting::kernel::create<TestType>();
    auto &ctx = kernel->context();
    ctx.particleTypes().add("A", 1.0);
    ctx.particleTypes().add("B", 1.0);
    ctx.topologyRegistry().addType("T");

    SECTION("Break due to large displacement with high rate, dimer") {
        // "A"-"A" with threshold 10 kBT
        // todo
        //auto &&breakingBonds = kernel->actions().breakBonds();
    }

    SECTION("Break due to large displacement with high rate, trimer") {

    }

    SECTION("Break due to large displacement with high rate, different types") {

    }

    SECTION("Break due to umbrella (box potential) pulling the bond apart (integration w/ forces and diffusion)") {

    }

}
