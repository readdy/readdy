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
 * << detailed description >>
 *
 * @file CPUKernelTest.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/actions/Actions.h>
#include <readdy/model/RandomProvider.h>

namespace {

TEST(CPUTestKernel, TestKernelLoad) {
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");

    kernel->context().boxSize() = {{10, 10, 10}};
    kernel->context().particleTypes().add("X", .55);
    kernel->context().periodicBoundaryConditions() = {{true, true, true}};
    kernel->context().reactions().addDecay("X decay", "X", .5);
    kernel->context().reactions().addFission("X fission", "X", "X", "X", .00, .5);

    auto &&integrator = kernel->actions().eulerBDIntegrator(1);
    auto &&neighborList = kernel->actions().updateNeighborList();
    auto &&forces = kernel->actions().calculateForces();
    auto &&reactions = kernel->actions().gillespie(1);

    auto pp_obs = kernel->observe().positions(1);
    auto connection = kernel->connectObservable(pp_obs.get());

    const int n_particles = 500;
    const auto typeId = kernel->context().particleTypes().idOf("X");
    std::vector<readdy::model::Particle> particlesToBeginWith{n_particles, {0, 0, 0, typeId}};
    readdy::log::debug("n_particles={}", particlesToBeginWith.size());
    kernel->stateModel().addParticles(particlesToBeginWith);

    neighborList->perform();
    for (size_t t = 0; t < 20; t++) {

        forces->perform();
        integrator->perform();
        neighborList->perform();

        reactions->perform();

        pp_obs->evaluate();
        readdy::log::debug("\tcurrently n particles: {}", pp_obs->getResult().size());
    }

    connection.disconnect();
}
}
