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
 * Catch sections (induced by SECTION, GIVEN, WHEN, THEN) cannot be used inside diverging program paths.
 * I.e. if rank 1 and rank 2 perform different branches of programs with matching barriers,
 * then a section inside that block would cause multiple execution runs.
 * This eventually results in non-matching barriers, undefined behaviors, and dead-locks.
 *
 * BAD: single instruction multiple data (SIMD) model is not valid anymore
 * if (rank == 0) {
 *     SECTION("BLUB") {
 *         REQUIRE(true);
 *     }
 * }
 *
 * GOOD:
 * SECTION("BLUB") {
 *     if (rank == 0) {
 *         REQUIRE(true);
 *     }
 * }
 *
 *
 * @file TestMain.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.05.19
 */


#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>

#include <readdy/plugin/KernelProvider.h>
#include <readdy/kernel/mpi/MPISession.h>
#include <readdy/plugin/Utils.h>


int main(int argc, char **argv) {
    readdy::kernel::mpi::MPISession mpisession(argc, argv);
    if (mpisession.rank() == 0) {
        readdy::log::console()->set_level(spdlog::level::info);
    } else {
        readdy::log::console()->set_level(spdlog::level::warn);
    }

    Catch::Session session;
    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) return returnCode;

    if (!session.config().listTestNamesOnly()) {
        const auto dir = readdy::plugin::utils::getPluginsDirectory();
        readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(dir);
    }

    return session.run();
}
