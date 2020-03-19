/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file TestPlugins.cpp
 * @author clonker
 * @date 07.03.16
 */

#include <catch2/catch.hpp>

#include <readdy/plugin/KernelProvider.h>
#include <readdy/testing/KernelMock.h>
#include <readdy/common/string.h>
#include <readdy/common/filesystem.h>

SCENARIO("Testing the kernel provider") {
    GIVEN("a kernel provider") {
        auto &provider = readdy::plugin::KernelProvider::getInstance();
        WHEN("trying to load a nonexisting problem") {
            THEN("an exception should be thrown") {
                REQUIRE_THROWS_AS(provider.create("foo"), std::invalid_argument);
            }
        }
        WHEN("loading an existing plugin") {
            auto kk_ptr = provider.create("SingleCPU");
            THEN("at least the names should match") {
                REQUIRE(kk_ptr->name() == "SingleCPU");
            }
        }
    }
}

TEST_CASE("Test Kernel ctor") {
    auto& provider = readdy::plugin::KernelProvider::getInstance();
    auto kernel = provider.create("SingleCPU");
    auto &ctx = kernel->context();
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    kernel->addParticle("A", readdy::Vec3(1, 0, 2));

    std::unordered_map<std::string, std::string> conversionsMap = {{"A", "B"}};
    //ctx.compartments().addSphere(conversionsMap, "kugelrund", readdy::Vec3(0, 0, 0), 10., false);
    ctx.topologyRegistry().addType("T");
    std::cout << ctx.topologyRegistry().isSpatialReactionType("T") << std::endl;
}
