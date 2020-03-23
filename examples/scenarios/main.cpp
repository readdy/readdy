/**
 * @file main.cpp
 * @brief Run scenarios for benchmarking, save output to json file
 * @author chrisfroe
 * @date 23.03.20
 */

#include <fstream>
#include <json.hpp>
#include "Scenarios.h"

using Json = nlohmann::json;
using Kernel = readdy::kernel::scpu::SCPUKernel;
//using Kernel = readdy::kernel::cpu::CPUKernel;
// MPI Kernel cannot be run here, it needs a MPISession and is called by mpirun,

int main(int argc, char **argv) {
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

    std::vector<std::unique_ptr<readdy::benchmark::Scenario<Kernel>>> scenarios;
    scenarios.push_back(std::make_unique<readdy::benchmark::FreeDiffusion<Kernel>>());

    for (const auto &s : scenarios) {
        Json json = s->run();

        json["scenarioName"] = s->name();
        json["scenarioDescription"] = s->description();
        json["kernel"];

        std::string filename{s->name()};
        std::string path = outDir + filename;

        if (not json.empty()) {
            std::ofstream stream(path, std::ofstream::out | std::ofstream::trunc);
            stream << json << std::endl;
        }
    }

    return 0;
}
