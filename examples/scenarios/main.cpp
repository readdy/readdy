/**
 * @file main.cpp
 * @brief Run scenarios for benchmarking or profiling, save output to json file
 * @author chrisfroe
 * @date 23.03.20
 */

#include <fstream>
#include <json.hpp>
#include <readdy/plugin/Utils.h>
#include <iomanip>
#include "Scenarios.h"

using Json = nlohmann::json;
namespace perf = readdy::performance;

/**
 * Arguments (all optional):
 * --outdir (where the output files will be stored)
 * --version (`git describe` when built from source or `conda list --json readdy` when using installed)
 * --cpu
 * --machine
 * --author
 * --prefix (user chosen, hinting at the purpose of this measurement. E.g. 'benchmark' or 'neighborlist-impl-3')
 *
 * Output files are written to outdir, each scenario to a
 * different file with name "${prefix}-${scenarioName}-${datetime}.json".
 * All given info (except outdir) will be embedded in the output file.
 *
 * Example for development:
 * run_readdy_scenarios \
 *  --outdir="/tmp" \
 *  --version="$(cd /path/to/readdy && git describe --always)" \
 *  --cpu="$(lscpu --json | tr -d '\n')" \
 *  --machine="$(hostname)" \
 *  --author="chris" \
 *  --prefix="playground"
 *
 *  Note: MPI Kernel cannot be run here, it needs a MPISession and is called by mpirun.
 */
int main(int argc, char **argv) {
    readdy::log::set_level(spdlog::level::info);

    // parse argument strings
    auto outdir = perf::getOption(argc, argv, "--outdir=", "/tmp/");
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

    auto pluginsDirectory = readdy::plugin::utils::getPluginsDirectory();
    readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(pluginsDirectory);

    // which scenarios shall be run
    std::vector<std::unique_ptr<readdy::performance::Scenario>> scenarios;
    {
        scenarios.push_back(std::make_unique<readdy::performance::FreeDiffusion>("SingleCPU", 100000));
        std::size_t load = readdy::readdy_default_n_threads();
        scenarios.push_back(std::make_unique<readdy::performance::DiffusionPairPotential>(
                "CPU", perf::WeakScalingGeometry::stick, load, 13.));
        scenarios.push_back(std::make_unique<readdy::performance::DiffusionPairPotential>(
                "SingleCPU", perf::WeakScalingGeometry::stick, load, 13.));
    }

    // run the scenarios, and write output
    for (const auto &s : scenarios) {
        readdy::log::info("Run scenario {} -- {}", s->name(), s->description());

        Json out;
        out["result"] = s->run();
        out["info"] = info;

        out["scenarioName"] = s->name();
        out["scenarioDescription"] = s->description();

        std::string filename = fmt::format("{}{}-{}-{}.json",
                prefix.empty() ? "" : prefix + "-", s->name(), time, perf::randomString());
        std::string path = outdir + filename;

        if (not out.empty()) {
            std::ofstream stream(path, std::ofstream::out | std::ofstream::trunc);
            stream << out << std::endl;
        }
    }

    return 0;
}
