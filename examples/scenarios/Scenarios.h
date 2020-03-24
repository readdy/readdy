/**
 * The goal here is to implement some reproducible scenarios e.g. for benchmark purpose.
 *
 * A note on readdy performance in reaction-diffusion problems:
 * readdy focuses on crowded cytosolic environments, diffusion-influenced kinetics and interactions via potentials.
 * Thus readdy should perform optimal for such cases. E.g. it is expected/obvious that an iPRD simulation with Brownian
 * Dynamics samples rather inefficient when simulating well-mixed kinetics and/or very dilute systems.
 *
 * Benchmark cases:
 * - Free diffusion
 * - Stationary distribution in an external potential -> diffusion in double well along one dimension
 * - Thermodynamics of LJ suspension -> diffusion subject to pair-interactions
 * - Michaelis-Menten kinetics -> irreversible reaction-diffusion without forces in the reaction-limit (well-mixed)
 * - A+B<-->C with LJ -> reversible reaction-diffusion with forces, can do diffusion-influenced
 *
 * @file Scenarios.h
 * @brief Generic kernel-independent scenarios to benchmark readdy performance
 * @author chrisfroe
 * @date 28.02.20
 */

#pragma once

#include <readdy/model/Kernel.h>
#include <readdy/api/Simulation.h>
#include <readdy/api/KernelConfiguration.h>
#include <readdy/common/Timer.h>
#include <utility>

using Json = nlohmann::json;
namespace rnd = readdy::model::rnd;

namespace readdy::performance {


class Scenario {
    std::string _description{};
    std::string _name{};
public:
    explicit Scenario(std::string name, std::string descr) : _name(std::move(name)), _description(std::move(descr)) {}

    const std::string &description() { return _description; }

    const std::string &name() { return _name; }

    virtual Json run() = 0;
};


class FreeDiffusion : public Scenario {
    std::string _kernelName;
    std::size_t _nSteps = 1000;
public:
    explicit FreeDiffusion(std::string kernelName) : Scenario(
            "FreeDiffusion"+kernelName,
            "Particles diffusing freely in periodic box without any interaction"),
            _kernelName(kernelName) {}

    static readdy::model::Context context() {
        readdy::model::Context ctx;
        ctx.boxSize() = {20., 20., 20.};
        ctx.particleTypes().add("A", 1.);
        return ctx;
    }

    Json simulate(std::size_t nParticles) {
        readdy::Simulation sim(_kernelName, context());
        auto &box = sim.context().boxSize();
        for (std::size_t i = 0; i < nParticles; ++i) {
            sim.addParticle("A",
                    rnd::uniform_real() * box[0] - 0.5 * box[0],
                    rnd::uniform_real() * box[1] - 0.5 * box[1],
                    rnd::uniform_real() * box[2] - 0.5 * box[2]);
        }

        {
            readdy::util::Timer t("totalSimulation");
            sim.run(_nSteps, 0.01);
        }

        auto perf = Json::parse(readdy::util::Timer::perfToJsonString());
        readdy::util::Timer::clear();
        return perf;
    }

    Json run() override {
        std::vector<std::size_t> numbers = {
                1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000
        };

        std::vector<Json> results;
        for (std::size_t iid = 0; iid < 10; ++iid) {
            for (auto n : numbers) {
                Json currentResult;
                currentResult["performance"] = simulate(n);
                currentResult["n"] = n;
                currentResult["iid"] = iid;
                results.push_back(currentResult);
            }
        }
        Json result = results;
        result["context"] = context().describe();
        result["nSteps"] = _nSteps;
        result["kernelName"] = _kernelName;
        result["readdy_default_n_threads"] = readdy_default_n_threads();

        return result;
    }
};

}
