/**
 * The goal here is _not_ to strictly measure scaling behavior but having fixed, reproducible benchmark cases.
 * ReaDDy focuses on crowded cytosolic environments, diffusion-influenced kinetics and interactions via potentials.
 * Thus ReaDDy should perform optimal for such cases. E.g. it is expected/obvious that an iPRD simulation performs
 * rather poorly when simulating well-mixed kinetics.
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

namespace readdy::benchmark {

using KernelPointer = plugin::KernelProvider::kernel_ptr;

template<typename Kernel>
class Scenario {
    std::string _description{};
    std::string _name{};
public:
    explicit Scenario(std::string name, std::string descr) : _name(std::move(name)), _description(std::move(descr)) {}

    const std::string &description() { return _description; }

    const std::string &name() { return _name; }

    virtual Json run() = 0;

    // todo if Kernel is a string, ask the KernelProvider to provide a kernel pointer
    // todo otherwise directly create
    KernelPointer createKernel() {
        return {std::make_unique<Kernel, readdy::plugin::KernelDeleter>()};
    }
};

template<typename Kernel>
class FreeDiffusion : public Scenario<Kernel> {
public:
    FreeDiffusion() : Scenario<Kernel>(
            "FreeDiffusion",
            "Particles diffusing freely in periodic box without any interaction") {}

    Json run() override {
        readdy::model::Context ctx;

        ctx.boxSize() = {20., 20., 20.};
        ctx.particleTypes().add("A", 1.);
        ctx.particleTypes().add("B", 1.);

        KernelPointer kernel(this->createKernel());
        readdy::Simulation sim(std::move(kernel), ctx);

        // todo sim

        Json result;
        result["benchmark"] = Json::parse(readdy::util::Timer::perfToJsonString());
        readdy::util::Timer::clear();
        return result;
    }
};

}
