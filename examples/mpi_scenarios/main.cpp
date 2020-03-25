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
 * @file main.cpp
 * @brief Run performance scenarios for the MPI kernel
 * @author chrisfroe
 * @date 28.05.19
 */

#include <readdy/kernel/mpi/MPISession.h>
#include <fstream>
#include <iomanip>
#include "MPIScenarios.h"

using Json = nlohmann::json;
namespace rkm = readdy::kernel::mpi;
namespace perf = readdy::performance;

int main(int argc, char **argv) {
    // MPI_Init will modify argc, argv such that they behave ''normal'' again, i.e. without the mpirun arguments
    rkm::MPISession mpiSession(argc, argv);

    readdy::log::set_level(spdlog::level::trace);

    // parse argument strings
    auto outdir = perf::getOption(argc, argv, "--outdir=", "/tmp");
    auto version = perf::getOption(argc, argv, "--version=", "no version info provided");
    auto cpuinfo = perf::getOption(argc, argv, "--cpu=", "no cpu info provided");
    auto machine = perf::getOption(argc, argv, "--machine=", "no machine name provided");
    auto author = perf::getOption(argc, argv, "--author=", "nobody");
    auto prefix = perf::getOption(argc, argv, "--prefix=", "");

    // necessary argument checking
    if (not(readdy::util::fs::exists(outdir) and readdy::util::fs::is_directory(outdir))) {
        throw std::invalid_argument(
                fmt::format("Target output directory {} does not exist or is no directory.", outdir));
    }

    // gather miscellaneous information
    auto time = perf::datetime();
    Json info;
    info["datetime"] = time;
    info["version"] = version;
    if (Json::accept(cpuinfo)) {
        Json cpuJson = Json::parse(cpuinfo);
        info["cpu"] = cpuJson;
    } else {
        info["cpu"] = cpuinfo;
    }
    info["machine"] = machine;
    info["author"] = author;
    info["prefix"] = prefix;

    info["rank"] = mpiSession.rank();
    info["worldSize"] = mpiSession.worldSize();
    info["processorName"] = mpiSession.processorName();

    // which scenarios shall be run
    std::vector<std::unique_ptr<perf::Scenario>> scenarios;
    //scenarios.push_back(std::make_unique<perf::MPIDistributeParticles>());
    scenarios.push_back(std::make_unique<perf::MPIDiffusionPairPotential>(perf::WeakScalingGeometry::stick));

    // run the scenarios, and write output
    for (const auto &s : scenarios) {
        readdy::log::info("rank={}, Run scenario {} -- {}", mpiSession.rank(), s->name(), s->description());
        Json out;
        rkm::MPISession::barrier();
        out["result"] = s->run();
        out["info"] = info;

        out["scenarioName"] = s->name();
        out["scenarioDescription"] = s->description();

        std::string filename = fmt::format("{}-{}-rank-{}-ws-{}.json", s->name(), time, mpiSession.rank(), mpiSession.worldSize());
        if (!prefix.empty()) {
            filename.insert(0, prefix + "-");
        }

        std::string path = outdir + filename;

//        if (not out.empty()) {
//            std::ofstream stream(path, std::ofstream::out | std::ofstream::trunc);
//            stream << out << std::endl;
//        }
    }

    readdy::log::info("rank={}, Done with scenarios", mpiSession.rank());
    return 0;
}
