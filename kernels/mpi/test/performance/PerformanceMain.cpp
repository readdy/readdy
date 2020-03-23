/********************************************************************
 * Copyright © 2020 Computational Molecular Biology Group,          *
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
 * Run performance scenarios (weak/strong scaling, different systems) for the MPI kernel.
 * Mainly measure runtime, thus no unit test framework required here.
 *
 * @file PerformanceMain.cpp
 * @brief Run performance scenarios for the MPI kernel
 * @author chrisfroe
 * @date 28.05.19
 */

#include <readdy/plugin/KernelProvider.h>
#include <readdy/kernel/mpi/MPISession.h>
#include <fstream>
#include "MPIScenarios.h"

using Json = nlohmann::json;
namespace rkm = readdy::kernel::mpi;

int main(int argc, char **argv) {
    // MPI_Init will modify argc, argv such that they behave ''normal'' again, i.e. without the mpirun arguments
    readdy::kernel::mpi::MPISession mpiSession(argc, argv);

    std::string outDir;
    if (argc > 1) {
        outDir = std::string(argv[1]);
    } else {
        readdy::log::warn("Using home directory as output");
        outDir = "~/";
    }

    std::vector<std::unique_ptr<readdy::benchmark::Scenario<rkm::MPIKernel>>> scenarios;
    scenarios.push_back(std::make_unique<rkm::benchmark::DistributeParticles>());

    for (const auto &s : scenarios) {
        mpiSession.barrier();
        auto json = s->run();

        json["rank"] = mpiSession.rank();
        json["worldSize"] = mpiSession.worldSize();
        json["processorName"] = mpiSession.processorName();
        json["scenarioName"] = s->name();
        json["scenarioDescription"] = s->description();

        std::string filename{s->name() + "_rank_" + std::to_string(mpiSession.rank())};
        std::string path = outDir + filename;

        if (not json.empty()) {
            std::ofstream stream(path, std::ofstream::out |std::ofstream::trunc);
            stream << json << std::endl;
        }
    }

    return 0;
}
