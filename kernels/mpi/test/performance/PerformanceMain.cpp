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
 * Run performance scenarios for the MPI kernel. As these do validate the correctness of readdy but mainly
 * measure runtime, this does not need the unit test library.
 *
 * @file PerformanceMain.cpp
 * @brief Run performance scenarios for the MPI kernel
 * @author chrisfroe
 * @date 28.05.19
 */

#include <readdy/plugin/KernelProvider.h>
#include <readdy/kernel/mpi/MPISession.h>
#include <readdy/plugin/Utils.h>
#include <fstream>
#include "Scenarios.h"

using json = nlohmann::json;
namespace perf = readdy::kernel::mpi::performance;

int main(int argc, char **argv) {
    // MPI_Init will modify argc, argv such that they behave ''normal'' again, i.e. without the mpirun arguments
    MPISession mpiSession(argc, argv);

    {
        const auto dir = readdy::plugin::utils::getPluginsDirectory();
        readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(dir);
    }

    std::string outDir;
    if (argc > 1) {
        outDir = std::string(argv[1]);
    } else {
        readdy::log::warn("Using home directory as output");
        outDir = "~/";
    }

    std::vector<std::unique_ptr<perf::Scenario>> scenarios;
    scenarios.push_back(std::make_unique<perf::DistributeParticles>());
    scenarios.push_back(std::make_unique<perf::DistributeParticles>());

    for (const auto &s : scenarios) {
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

        MPI_Barrier(MPI_COMM_WORLD);
    }

    return 0;
}
