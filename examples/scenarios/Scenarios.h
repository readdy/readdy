/**
 * The goal here is to implement some reproducible scenarios for performance measurement
 * Not necessarily benchmark, more profiling, bottleneck identification, AB comparisons.
 *
 * A note on readdy performance in reaction-diffusion problems:
 * readdy focuses on crowded cytosolic environments, diffusion-influenced kinetics and interactions via potentials.
 * Thus readdy should perform optimal for such cases. E.g. it is expected/obvious that an iPRD simulation with Brownian
 * Dynamics samples rather inefficient when simulating well-mixed kinetics and/or very dilute systems.
 *
 * @file Scenarios.h
 * @brief Generic kernel-independent scenarios to measure readdy performance
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

std::string randomString(std::size_t n = 6) {
    assert(n > 0);
    assert(n < 32);
    std::string result("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
    std::random_device rd;
    std::mt19937 generator(rd());
    std::shuffle(result.begin(), result.end(), generator);
    return result.substr(0, n);
}

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
    std::size_t _nParticles = 1000;
public:
    explicit FreeDiffusion(const std::string& kernelName, std::size_t nParticles) : Scenario(
            "FreeDiffusion"+kernelName,
            "Particles diffusing freely in periodic box without any interaction"),
            _kernelName(kernelName), _nParticles(nParticles) {}

    Json run() override {
        readdy::model::Context ctx;
        ctx.boxSize() = {20., 20., 20.};
        ctx.particleTypes().add("A", 1.);
        readdy::Simulation sim(_kernelName, ctx);
        auto &box = sim.context().boxSize();
        for (std::size_t i = 0; i < _nParticles; ++i) {
            sim.addParticle("A",
                            rnd::uniform_real() * box[0] - 0.5 * box[0],
                            rnd::uniform_real() * box[1] - 0.5 * box[1],
                            rnd::uniform_real() * box[2] - 0.5 * box[2]);
        }

        readdy::util::Timer::clear();
        {
            readdy::util::Timer t("totalSimulation");
            sim.run(_nSteps, 0.01);
        }

        Json result;
        result["context"] = ctx.describe();
        result["nSteps"] = _nSteps;
        result["nParticles"] = _nParticles;
        result["kernelName"] = _kernelName;
        result["readdy_default_n_threads"] = readdy_default_n_threads();
        result["performance"] = Json::parse(readdy::util::Timer::perfToJsonString());
        readdy::util::Timer::clear();
        return result;
    }
};

/** Scale the box based on given load in different ways */
enum WeakScalingGeometry {
    stick, // scale box only in x.
    slab, // scale box in x and y.
    cube // scale box in x, y and z.
};

inline std::ostream &operator<<(std::ostream& os, WeakScalingGeometry mode) {
    switch(mode) {
        case stick: os << "stick"; break;
        case slab: os << "slab"; break;
        case cube: os << "cube"; break;
    }
    return os;
}

class DiffusionPairPotential : public Scenario {
    std::string _kernelName;
    WeakScalingGeometry _mode;
    // scales the amount of particles as well as the volume
    std::size_t _nLoad;
    // hyperparameter determining the volume and thus the final number of particles
    // describes the ratio of box length and interaction distance
    readdy::scalar _edgeLengthOverInteractionDistance;
public:
    DiffusionPairPotential(std::string kernelName, WeakScalingGeometry mode, std::size_t nLoad = 1, readdy::scalar edgeLengthOverInteractionDistance = 5.) : Scenario(
        "DiffusionPairPotential"+kernelName,
        "Diffusion of particles with constant density, scale box according to mode"),
        _mode(mode), _kernelName(kernelName), _nLoad(nLoad), _edgeLengthOverInteractionDistance(edgeLengthOverInteractionDistance) {
        assert(nLoad > 0);
        assert(edgeLengthOverInteractionDistance > 0.);
    }

    Json run() override {
        std::size_t nParticles;
        std::size_t nParticlesPerLoad;
        std::array<readdy::scalar, 3> box{};
        readdy::scalar interactionDistance = 2.;
        // radius that is used to calculate the volume occupation, if the particles were hard-spheres
        readdy::scalar assumedRadius = interactionDistance / 2.;
        if (_mode == stick) {
            // will lead to 60% volume occupation when one particle has an associated radius of halo/2.
            nParticlesPerLoad = 0.6 * std::pow(_edgeLengthOverInteractionDistance, 3)
                                            / (4./3. * readdy::util::numeric::pi<scalar>() * std::pow(assumedRadius / interactionDistance, 3));
            nParticles = _nLoad * nParticlesPerLoad;
            readdy::scalar boxLength = _edgeLengthOverInteractionDistance * interactionDistance;
            box = {_nLoad * boxLength, boxLength, boxLength};
        } else if (_mode == cube) {
            nParticlesPerLoad = 0.6 * std::pow(_edgeLengthOverInteractionDistance, 3)
                                / (4./3. * readdy::util::numeric::pi<scalar>() * std::pow(assumedRadius/interactionDistance,3));
            nParticles = _nLoad * nParticlesPerLoad;
            scalar boxLength = _edgeLengthOverInteractionDistance * interactionDistance;
            scalar factor = std::cbrt(static_cast<scalar>(_nLoad));
            box = {factor * boxLength, factor * boxLength, factor * boxLength};
        } else {
            throw std::invalid_argument(fmt::format("Unknown scaling mode {}", _mode));
        }
        readdy::log::info("nParticlesPerLoad {}, nParticles {}", nParticlesPerLoad, nParticles);

        readdy::model::Context ctx;
        ctx.boxSize() = box;
        ctx.particleTypes().add("A", 1.);
        ctx.potentials().addHarmonicRepulsion("A", "A", 10., interactionDistance);

        auto kernel = readdy::plugin::KernelProvider::getInstance().create(_kernelName);
        kernel->context() = ctx;

        auto idA = kernel->context().particleTypes().idOf("A");
        std::vector<readdy::model::Particle> particles;
        for (std::size_t i = 0; i < nParticles; ++i) {
            auto x = readdy::model::rnd::uniform_real() * ctx.boxSize()[0] - 0.5 * ctx.boxSize()[0];
            auto y = readdy::model::rnd::uniform_real() * ctx.boxSize()[1] - 0.5 * ctx.boxSize()[1];
            auto z = readdy::model::rnd::uniform_real() * ctx.boxSize()[2] - 0.5 * ctx.boxSize()[2];
            particles.emplace_back(x, y, z, idA);
        }

        auto addParticles = kernel->actions().addParticles(particles);
        {
            readdy::util::Timer t("addParticles");
            addParticles->perform();
        }

        readdy::scalar timeStep = 0.01;
        std::size_t nSteps = 1000;
        auto integrator = kernel->actions().eulerBDIntegrator(timeStep);
        auto forces = kernel->actions().calculateForces();
        auto createNL = kernel->actions().createNeighborList(kernel->context().calculateMaxCutoff());
        auto neighborList = kernel->actions().updateNeighborList();

        createNL->perform();
        neighborList->perform();
        forces->perform();
        for (size_t t = 1; t < nSteps + 1; t++) {
            readdy::util::Timer tStep("complete timestep");
            {
                readdy::util::Timer t1("integrator");
                integrator->perform();
            }
            {
                readdy::util::Timer t2("neighborList");
                neighborList->perform();
            }
            {
                readdy::util::Timer t3("forces");
                forces->perform();
            }
        }

        std::size_t nProcessors = _kernelName == "SingleCPU" ? 1 : readdy_default_n_threads();
        Json result;
        result["context"] = kernel->context().describe();
        result["performance"] = Json::parse(readdy::util::Timer::perfToJsonString());
        result["nLoad"] = _nLoad;
        result["nParticles"] = nParticles;
        result["nParticlesPerLoad"] = nParticlesPerLoad;
        result["nProcessors"] = nProcessors;
        result["nParticlesPerProcessor"] = static_cast<scalar>(nParticles) / static_cast<scalar>(nProcessors);
        result["kernelName"] = _kernelName;
        result["mode"] = fmt::format("{}", _mode);
        result["edgeLengthOverInteractionDistance"] = _edgeLengthOverInteractionDistance;
        readdy::util::Timer::clear();
        return result;
    }
};

}
