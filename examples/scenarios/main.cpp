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

std::string datetime() {
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::time_t nowTime = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::gmtime(&nowTime), "%F-%T");
    return ss.str();
}

std::string getOption(int argc, char **argv, const std::string &option, const std::string &defaultValue = "") {
    std::string value;
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find(option) == 0) { // C++20 has starts_with
            value = arg.substr(option.size());
            return value;
        }
    }
    return defaultValue;
}

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
 *  --version="$(cd /home/chris/workspace/readdy && git describe --always)" \
 *  --cpu="$(lscpu --json | tr -d '\n')" \
 *  --machine="$(hostname)" \
 *  --author="chris" \
 *  --prefix="playground"
 *
 *  Note: MPI Kernel cannot be run here, it needs a MPISession and is called by mpirun.
 */
int main(int argc, char **argv) {
    // parse argument strings
    auto outdir = getOption(argc, argv, "--outdir=", "/tmp");
    auto version = getOption(argc, argv, "--version=", "no version info provided");
    auto cpuinfo = getOption(argc, argv, "--cpu=", "no cpu info provided");
    auto machine = getOption(argc, argv, "--machine=", "no machine name provided");
    auto author = getOption(argc, argv, "--author=", "nobody");
    auto prefix = getOption(argc, argv, "--prefix=", "");

    // necessary argument checking
    if (not(readdy::util::fs::exists(outdir) and readdy::util::fs::is_directory(outdir))) {
        throw std::invalid_argument(
                fmt::format("Target output directory {} does not exist or is no directory.", outdir));
    }

    // gather miscellaneous information
    auto time = datetime();
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

    auto pluginsdir = readdy::plugin::utils::getPluginsDirectory();
    readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(pluginsdir);

    // which scenarios shall be run
    std::vector<std::unique_ptr<readdy::performance::Scenario>> scenarios;
    scenarios.push_back(std::make_unique<readdy::performance::FreeDiffusion>("SingleCPU"));
    scenarios.push_back(std::make_unique<readdy::performance::FreeDiffusion>("CPU"));

    // run the scenarios, and write output
    for (const auto &s : scenarios) {
        Json out;
        out["result"] = s->run();
        out["info"] = info;

        out["scenarioName"] = s->name();
        out["scenarioDescription"] = s->description();

        std::string filename = s->name() + "-" + time;
        if (!prefix.empty()) {
            filename.insert(0, prefix + "-");
        }

        std::string path = outdir + filename + ".json";

        if (not out.empty()) {
            std::ofstream stream(path, std::ofstream::out | std::ofstream::trunc);
            stream << out << std::endl;
        }
    }

    return 0;
}
