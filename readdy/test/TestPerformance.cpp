/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * These tests should cover different scenarios of system configurations. They should provoke very different
 * computational efforts.
 * We define the following unitless quantities (N,L,S,R), which will help to compare different scenarios:
 *  - N ~ number of particles. For reactive systems this is the initial number of particles. Since we do measure for a few timesteps
 *    only, this is close to the average number of particles.
 *  - L = V / h**3 ~ system size, i.e. fraction of system volume and the largest relevant lengthscale. For a cubic system V**1/3 would
 *    be the edge length of the box. With repulsing particles of radius r, the largest relevant lengthscale
 *    would be h=2r. V is the volume _accessible_ to particles. Note that h is the largest distance of two particles
 *    when an interaction has to be considered. Thus h can also be a reaction distance.
 *  - S = sqrt(2*D*Delta t) / h ~ relative displacement of particles, fraction of distance traveled within
 *    one timestep (of width Delta t) divided by the largest relevant lengthscale. D is the largest diffusion constant.
 *    Usually S should be smaller than 1.
 *  - R = k*Delta t ~ reactivity. When k is the largest reaction rate, R should usually be smaller than 1.
 *
 *  - Derived from N and L can be the dimensionless density rho = N / L. When there is only a single relevant lengthscale h
 *    or all lengthscales are narrowly distributed around h, the following statements apply:
 *      - when rho = 1, the system is rather dense (every particle is assigned with a cube of volume h**3)
 *      - when rho ~= 2, the system is closely packed (every particles is a sphere of radius h/2)
 *      - when rho < 1, the system is dilute
 *      - when rho >> 1, either lengthscales are widely distributed (then rho is ambiguous) or the system is ill-parametrized
 *
 * There are different combinations of scenarios where our algorithms can scale very differently:
 *  - Purely reactive system vs. purely collisive system vs. both reactive and collisive
 *  - Uniformly distributed vs. spatially separated
 *  - Homogeneously sized particles vs. heterogene mixture
 *
 * If you consider ONE system (fixed N and L) and want to compare different simulation algorithms with
 * varying timestep Delta t, then the simulation efficiency is of central interest:
 *  - eta = simulated time / computation time ~ simulation efficieny
 *
 * For example increasing the timestep would increase eta but R and S will also increase. If R or S come into range
 * of ~1, this hints towards a sampling problem.
 *
 * @file TestPerformance.cpp
 * @brief Performance test with different scenarios. Per default not included in tests.
 * @author clonker
 * @author chrisfroe
 * @date 11.07.16
 */

#include <gtest/gtest.h>

#include <readdy/testing/Utils.h>
#include <readdy/common/Timer.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/Utils.h>
#include <readdy/io/File.h>
#include <readdy/api/Simulation.h>

namespace {

using file_t = readdy::io::File;
using update_neighbor_list_t = readdy::model::actions::UpdateNeighborList;

const std::string SKIN_FACTOR = "skin";
const std::string NUMBERS_FACTOR = "numbers";
const std::string BOXLENGTHS_FACTOR = "boxlengths";

/**
 * A performance scenario is an instance of a simulation that is run once for a few timesteps to record
 * the duration of force calculation, diffusion integration, neighborlist updating and performing reactions.
 */
class PerformanceScenario {
public:
    PerformanceScenario(const std::string &kernelName, readdy::scalar timeStep, readdy::scalar skin) :
            kernel(readdy::plugin::KernelProvider::getInstance().create(kernelName)), timeStep(timeStep),
            clearNeighborList(kernel->createAction<update_neighbor_list_t>(update_neighbor_list_t::Operation::clear)) {
        neighborList = kernel->createAction<update_neighbor_list_t>(update_neighbor_list_t::Operation::create, skin);
    }

    template<typename ReactionScheduler=readdy::model::actions::reactions::Gillespie>
    void perform(const readdy::time_step_type steps = 5, bool verbose = false) {
        this->configure();
        const auto result = runPerformanceTest<ReactionScheduler>(steps, verbose);
        timeForces = std::get<0>(result);
        timeIntegrator = std::get<1>(result);
        timeNeighborList = std::get<2>(result);
        timeReactions = std::get<3>(result);
        hasPerformed = true;
    };

    readdy::scalar getTimeForces() const {
        if (hasPerformed) return timeForces;
        else throw std::runtime_error("scenario has not performed yet");
    }

    readdy::scalar getTimeIntegrator() const {
        if (hasPerformed) return timeIntegrator;
        else throw std::runtime_error("scenario has not performed yet");
    }

    readdy::scalar getTimeNeighborlist() const {
        if (hasPerformed) return timeNeighborList;
        else throw std::runtime_error("scenario has not performed yet");
    }

    readdy::scalar getTimeReactions() const {
        if (hasPerformed) return timeReactions;
        else throw std::runtime_error("scenario has not performed yet");
    }

    readdy::scalar getParticleNumber() const { return particleNumber; }

    readdy::scalar getSystemSize() const { return systemSize; }

    readdy::scalar getRelativeDisplacement() const { return relativeDisplacement; }

    readdy::scalar getReactivity() const { return reactivity; }

protected:
    /**
     * Run a simulation for a given number of timesteps. All configuration is included in the kernel.
     * Which reaction scheduler is used is determined via template parameters
     *
     * @param steps number of timesteps that will be performed
     * @param verbose decides if to print information on context and performance
     * @return computation time per timestep for {forces, integrator, neighborlist, reactions}
     */
    template<typename ReactionScheduler=readdy::model::actions::reactions::Gillespie>
    std::array<readdy::scalar, 4>
    runPerformanceTest(readdy::time_step_type steps, const bool verbose = false) {
        auto &&integrator = kernel->createAction<readdy::model::actions::EulerBDIntegrator>(timeStep);
        auto &&forces = kernel->createAction<readdy::model::actions::CalculateForces>();
        auto &&reactions = kernel->createAction<ReactionScheduler>(timeStep);
        kernel->getKernelContext().configure(verbose);

        neighborList->perform();
        {
            forces->perform(true, "forces");
        }
        for (readdy::time_step_type t = 0; t < steps; ++t) {
            if (verbose) {
                readdy::log::debug("----------");
                readdy::log::debug("t = {}", t);
            }
            integrator->perform(true, "integrator");
            neighborList->perform(true, "neighbor list 1");
            reactions->perform(true, "reactions");
            neighborList->perform(true, "neighbor list 2");
            forces->perform(true, "forces");
        }
        clearNeighborList->perform();
        readdy::log::critical("DONE!");
        using timer = readdy::util::Timer;
        timeForces = timer::times().at("forces");
        timeIntegrator = timer::times().at("integrator");
        timeNeighborList = timer::times().at("neighbor list 1") + timer::times().at("neighbor list 2");
        timeReactions = timer::times().at("reactions");
        timeForces /= timer::counts().at("forces");
        timeIntegrator /= timer::counts().at("integrator");
        timeNeighborList /= timer::counts().at("neighbor list 1") + timer::counts().at("neighbor list 2");
        timeReactions /= timer::counts().at("reactions");
        if (verbose) {
            std::cout << "--------------------------------------------------------------" << std::endl;
            std::cout << "Average time for calculating forces: " << timeForces << std::endl;
            std::cout << "Average time for the integrator:     " << timeIntegrator << std::endl;
            std::cout << "Average time for the neighbor list:  " << timeNeighborList << std::endl;
            std::cout << "Average time for handling reactions: " << timeReactions << std::endl;
            std::cout << "--------------------------------------------------------------" << std::endl;
        }
        const std::array<readdy::scalar, 4> result = {timeForces, timeIntegrator, timeNeighborList, timeReactions};
        return result;

    }

    /** This is where the actual system gets set up: define reactions/potentials, add particles ... */
    virtual void configure() = 0;

    readdy::scalar timeForces = 0, timeIntegrator = 0, timeNeighborList = 0, timeReactions = 0; // main result
    readdy::scalar particleNumber = 0, systemSize = 0, relativeDisplacement = 0, reactivity = 0;
    readdy::scalar timeStep;

protected:
    // the (N,L,S,R) observables defined above
    readdy::plugin::KernelProvider::kernel_ptr kernel;
    std::unique_ptr<readdy::model::actions::UpdateNeighborList> neighborList;
    std::unique_ptr<readdy::model::actions::UpdateNeighborList> clearNeighborList;
    bool hasPerformed = false;
};

/**
 * The following scenario shall simulate a purely reactive system. Particles are uniformly distributed and
 * particle-radii/reaction-radii are homogeneous, i.e. there is only one length scale involved, the reaction radius = 4.5.
 * Reactions are A + B <--> C, with rates 'on' for the fusion and 'off' for the fission.
 * Parameter values resemble a biological cell in the cytosol where particles are proteins. To obtain real units the following
 * rescaling units are used:
 *  - energy is measured in 2.437 kJ/mol (which is KBT at room temperature)
 *  - lengths are measured in nm
 *  - time is measured in ns
 */
class ReactiveUniformHomogeneous : public PerformanceScenario {
public:
    static const std::string name;

    ReactiveUniformHomogeneous(const std::string &kernelName, const std::map<std::string, readdy::scalar> &factors)
            : PerformanceScenario(kernelName, .1,
                                  factors.find(SKIN_FACTOR) != factors.end() ? 4.5 * factors.at(SKIN_FACTOR) : -1) {
        numberA = static_cast<unsigned long>(numberA * factors.at(NUMBERS_FACTOR));
        numberC = static_cast<unsigned long>(numberC * factors.at(NUMBERS_FACTOR));
        boxLength *= factors.at(BOXLENGTHS_FACTOR);
        // set (N,L,S,R) observables
        particleNumber = 2 * numberA + numberC;
        systemSize = std::pow(boxLength, 3.) / std::pow(4.5, 3.);
        relativeDisplacement = std::sqrt(2 * diffusionConstants["A"] * timeStep) / 4.5;
        reactivity = rateOn * timeStep;
    };

    virtual void configure() override {
        auto &ctx = kernel->getKernelContext();
        ctx.setKBT(1.);
        ctx.setBoxSize(boxLength, boxLength, boxLength);
        ctx.setPeriodicBoundary(true, true, true);
        for (const std::string &type : {"A", "B", "C"}) {
            ctx.particle_types().add(type, diffusionConstants[type], radii[type]);
        }

        kernel->registerReaction<readdy::model::reactions::Fusion>("A+B->C", "A", "B", "C", rateOn, 4.5);
        kernel->registerReaction<readdy::model::reactions::Fission>("C->A+B", "C", "A", "B", rateOff, 4.5);

        /** Distribute particles uniformly in box */
        std::vector<readdy::model::Particle> particles;
        {
            particles.reserve(2 * numberA + numberC);
            auto uniform = [this]() {
                return readdy::model::rnd::uniform_real<readdy::scalar, std::mt19937>(-0.5 * boxLength, 0.5 * boxLength);
            };
            const auto &typeMapping = ctx.particle_types().type_mapping();
            const auto typeA = typeMapping.at("A");
            const auto typeB = typeMapping.at("B");
            const auto typeC = typeMapping.at("C");
            for (auto i = 0; i < numberA; ++i) {
                readdy::model::Particle particleA(uniform(), uniform(), uniform(), typeA);
                particles.push_back(particleA);
                readdy::model::Particle particleB(uniform(), uniform(), uniform(), typeB);
                particles.push_back(particleB);
            }
            for (auto i = 0; i < numberC; ++i) {
                readdy::model::Particle particleC(uniform(), uniform(), uniform(), typeC);
                particles.push_back(particleC);
            }
        }
        kernel->getKernelStateModel().addParticles(particles);
    }

private:
    unsigned long numberA = 500, numberC = 1800;
    readdy::scalar boxLength = 100., rateOn = 1e-3, rateOff = 5e-5;
    std::map<std::string, readdy::scalar> radii{{"A", 1.5},
                                                {"B", 3},
                                                {"C", 3.12}};
    std::map<std::string, readdy::scalar> diffusionConstants{{"A", .14},
                                                             {"B", .07},
                                                             {"C", .06}};
};

const std::string ReactiveUniformHomogeneous::name = "ReactiveUniformHomogeneous";

/**
 * The following scenario is the same as ReactiveUniformHomogeneous but adds repulsive forces
 */
class ReactiveCollisiveUniformHomogeneous : public PerformanceScenario {
public:
    static const std::string name;

    ReactiveCollisiveUniformHomogeneous(const std::string &kernelName, const std::map<std::string, readdy::scalar> &factors)
            : PerformanceScenario(kernelName, .1,
                                  factors.find(SKIN_FACTOR) != factors.end() ? 4.5 * factors.at(SKIN_FACTOR) : -1) {
        numberA = static_cast<unsigned long>(numberA * factors.at(NUMBERS_FACTOR));
        numberC = static_cast<unsigned long>(numberC * factors.at(NUMBERS_FACTOR));
        boxLength *= factors.at(BOXLENGTHS_FACTOR);
        // set (N,L,S,R) observables
        particleNumber = 2 * numberA + numberC;
        systemSize = std::pow(boxLength, 3.) / std::pow(4.5, 3.);
        relativeDisplacement = std::sqrt(2 * diffusionConstants["A"] * timeStep) / 4.5;
        reactivity = rateOn * timeStep;
    };

    virtual void configure() override {
        auto &ctx = kernel->getKernelContext();
        ctx.setKBT(1.);
        ctx.setBoxSize(boxLength, boxLength, boxLength);
        ctx.setPeriodicBoundary(true, true, true);
        for (const std::string &type : {"A", "B", "C"}) {
            ctx.particle_types().add(type, diffusionConstants[type], radii[type]);
        }
        kernel->registerReaction<readdy::model::reactions::Fusion>("A+B->C", "A", "B", "C", rateOn, 4.5);
        kernel->registerReaction<readdy::model::reactions::Fission>("C->A+B", "C", "A", "B", rateOff, 4.5);

        std::vector<std::pair<std::string, std::string>> pairs = {{"A", "B"},
                                                                  {"B", "C"},
                                                                  {"A", "C"}};
        for (const auto &typePair : pairs) {
            kernel->registerPotential<readdy::model::potentials::HarmonicRepulsion>(std::get<0>(typePair),
                                                                                    std::get<1>(typePair),
                                                                                    forceConstant);
        }

        /** Distribute particles uniformly in box */
        std::vector<readdy::model::Particle> particles;
        {
            particles.reserve(2 * numberA + numberC);
            auto uniform = [this]() {
                return readdy::model::rnd::uniform_real<readdy::scalar, std::mt19937>(-0.5 * boxLength, 0.5 * boxLength);
            };
            const auto &typeMapping = ctx.particle_types().type_mapping();
            const auto typeA = typeMapping.at("A");
            const auto typeB = typeMapping.at("B");
            const auto typeC = typeMapping.at("C");
            for (auto i = 0; i < numberA; ++i) {
                readdy::model::Particle particleA(uniform(), uniform(), uniform(), typeA);
                particles.push_back(particleA);
                readdy::model::Particle particleB(uniform(), uniform(), uniform(), typeB);
                particles.push_back(particleB);
            }
            for (auto i = 0; i < numberC; ++i) {
                readdy::model::Particle particleC(uniform(), uniform(), uniform(), typeC);
                particles.push_back(particleC);
            }
        }
        kernel->getKernelStateModel().addParticles(particles);
    }

private:
    unsigned long numberA = 500, numberC = 1800;
    readdy::scalar boxLength = 100., rateOn = 1e-3, rateOff = 5e-5, forceConstant = 10.;
    std::map<std::string, readdy::scalar> radii{{"A", 1.5},
                                                {"B", 3},
                                                {"C", 3.12}};
    std::map<std::string, readdy::scalar> diffusionConstants{{"A", .14},
                                                             {"B", .07},
                                                             {"C", .06}};
};

const std::string ReactiveCollisiveUniformHomogeneous::name = "ReactiveCollisiveUniformHomogeneous";

class CollisiveUniformHomogeneous : public PerformanceScenario {
public:
    static const std::string name;

    CollisiveUniformHomogeneous(const std::string &kernelName, const std::map<std::string, readdy::scalar> &factors)
            : PerformanceScenario(kernelName, .1,
                                  factors.find(SKIN_FACTOR) != factors.end() ? 4.5 * factors.at(SKIN_FACTOR) : -1) {
        numberA = static_cast<unsigned long>(numberA * factors.at(NUMBERS_FACTOR));
        numberC = static_cast<unsigned long>(numberC * factors.at(NUMBERS_FACTOR));
        boxLength *= factors.at(BOXLENGTHS_FACTOR);
        // set (N,L,S,R) observables
        particleNumber = 2 * numberA + numberC;
        systemSize = std::pow(boxLength, 3.) / std::pow(4.5, 3.);
        relativeDisplacement = std::sqrt(2 * diffusionConstants["A"] * timeStep) / 4.5;
        reactivity = 0;
    };

    virtual void configure() override {
        auto &ctx = kernel->getKernelContext();
        ctx.setKBT(1.);
        ctx.setBoxSize(boxLength, boxLength, boxLength);
        ctx.setPeriodicBoundary(true, true, true);
        for (const std::string &type : {"A", "B", "C"}) {
            ctx.particle_types().add(type, diffusionConstants[type], radii[type]);
        }
        std::vector<std::pair<std::string, std::string>> pairs = {{"A", "B"},
                                                                  {"B", "C"},
                                                                  {"A", "C"}};
        for (const auto &typePair : pairs) {
            kernel->registerPotential<readdy::model::potentials::HarmonicRepulsion>(std::get<0>(typePair),
                                                                                    std::get<1>(typePair),
                                                                                    forceConstant);
        }

        /** Distribute particles uniformly in box */
        std::vector<readdy::model::Particle> particles;
        {
            particles.reserve(2 * numberA + numberC);
            auto uniform = [this]() {
                return readdy::model::rnd::uniform_real<readdy::scalar, std::mt19937>(-0.5 * boxLength, 0.5 * boxLength);
            };
            const auto &typeMapping = ctx.particle_types().type_mapping();
            const auto typeA = typeMapping.at("A");
            const auto typeB = typeMapping.at("B");
            const auto typeC = typeMapping.at("C");
            for (auto i = 0; i < numberA; ++i) {
                readdy::model::Particle particleA(uniform(), uniform(), uniform(), typeA);
                particles.push_back(particleA);
                readdy::model::Particle particleB(uniform(), uniform(), uniform(), typeB);
                particles.push_back(particleB);
            }
            for (auto i = 0; i < numberC; ++i) {
                readdy::model::Particle particleC(uniform(), uniform(), uniform(), typeC);
                particles.push_back(particleC);
            }
        }
        kernel->getKernelStateModel().addParticles(particles);
    }

private:
    unsigned long numberA = 500, numberC = 1800;
    readdy::scalar boxLength = 100., forceConstant = 10.;
    std::map<std::string, readdy::scalar> radii{{"A", 1.5},
                                                {"B", 3},
                                                {"C", 3.12}};
    std::map<std::string, readdy::scalar> diffusionConstants{{"A", .14},
                                                             {"B", .07},
                                                             {"C", .06}};
};

const std::string CollisiveUniformHomogeneous::name = "CollisiveUniformHomogeneous";

template<typename scenario_t, typename ReactionScheduler=readdy::model::actions::reactions::Gillespie>
void scaleNumbersAndBoxsize(const std::string &kernelName) {
    /** Base values will be multiplied by factors. numbers[i] and boxlength[i] factors for same i will conserve particle density */
    const std::vector<readdy::scalar> numbers = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5,
                                                 3., 4., 5., 6., 7.};
    std::vector<readdy::scalar> boxLengths(numbers.size());
    {
        auto x = 0;
        std::generate(boxLengths.begin(), boxLengths.end(), [&numbers, &x] { return std::cbrt(numbers[x++]); });
    }
    const auto numbersSize = numbers.size();
    const auto boxLengthsSize = boxLengths.size();

    // results are performance times and the (N,L,S,R) observables
    readdy::scalar timeForces[numbersSize][boxLengthsSize];
    readdy::scalar timeIntegrator[numbersSize][boxLengthsSize];
    readdy::scalar timeNeighborList[numbersSize][boxLengthsSize];
    readdy::scalar timeReactions[numbersSize][boxLengthsSize];
    readdy::scalar particleNumber[numbersSize][boxLengthsSize];
    readdy::scalar systemSize[numbersSize][boxLengthsSize];
    readdy::scalar relativeDisplacement[numbersSize][boxLengthsSize];
    readdy::scalar reactivity[numbersSize][boxLengthsSize];

    for (auto n = 0; n < numbers.size(); ++n) {
        for (auto l = 0; l < boxLengths.size(); ++l) {
            std::map<std::string, readdy::scalar> factors;
            factors.emplace(std::make_pair(NUMBERS_FACTOR, numbers[n]));
            factors.emplace(std::make_pair(BOXLENGTHS_FACTOR, boxLengths[l]));
            scenario_t scenario(kernelName, factors);
            scenario.template perform<ReactionScheduler>(5);
            timeForces[n][l] = scenario.getTimeForces();
            timeIntegrator[n][l] = scenario.getTimeIntegrator();
            timeNeighborList[n][l] = scenario.getTimeNeighborlist();
            timeReactions[n][l] = scenario.getTimeReactions();
            particleNumber[n][l] = scenario.getParticleNumber();
            systemSize[n][l] = scenario.getSystemSize();
            relativeDisplacement[n][l] = scenario.getRelativeDisplacement();
            reactivity[n][l] = scenario.getReactivity();
        }
    }

    /** Write the result vectors to one dataset each, all in a single file */
    const std::string filename = "numbersAndBoxlengths_" + scenario_t::name + "_" + kernelName + ".h5";

    {
        file_t file(filename, file_t::Action::CREATE, file_t::Flag::OVERWRITE);
        auto inputGroup = file.createGroup("/input");
        inputGroup.write(NUMBERS_FACTOR, numbers);
        inputGroup.write(BOXLENGTHS_FACTOR, boxLengths);

        auto outputGroup = file.createGroup("/output");
        outputGroup.write("time_forces", {numbersSize, boxLengthsSize}, &timeForces[0][0]);
        outputGroup.write("time_integrator", {numbersSize, boxLengthsSize}, &timeIntegrator[0][0]);
        outputGroup.write("time_neighborlist", {numbersSize, boxLengthsSize}, &timeNeighborList[0][0]);
        outputGroup.write("time_reactions", {numbersSize, boxLengthsSize}, &timeReactions[0][0]);
        outputGroup.write("particle_number", {numbersSize, boxLengthsSize}, &particleNumber[0][0]);
        outputGroup.write("system_size", {numbersSize, boxLengthsSize}, &systemSize[0][0]);
        outputGroup.write("relative_displacement", {numbersSize, boxLengthsSize}, &relativeDisplacement[0][0]);
        outputGroup.write("reactivity", {numbersSize, boxLengthsSize}, &reactivity[0][0]);
    }
}

template<typename Scenario_t, typename ReactionScheduler=readdy::model::actions::reactions::Gillespie>
void scaleNumbersAndSkin(const std::string &kernelName, bool reducedNumbers) {
    readdy::log::debug("using reduced numbers: {}", reducedNumbers);
    /** Base values will be multiplied by factors. numbers[i] and boxlength[i] factors for same i will conserve particle density */
    std::vector<readdy::scalar> numbers;
    if (reducedNumbers) {
        numbers = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5,
                   3., 4., 5., 6., 7., 8., 9.,
                   10., 12., 15., 20., 40.};
    } else {
        numbers = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5,
                   3., 4., 5., 6., 7., 8., 9.,
                   10., 12., 15., 20., 40., 60, 80, 100};
    }
    std::vector<readdy::scalar> skinSizes;
    if (readdy::plugin::KernelProvider::getInstance().create(
            kernelName)->createAction<readdy::model::actions::UpdateNeighborList>()->supportsSkin()) {
        skinSizes = {.0, .1, .2, .3, .4, .5, 1, 2, 3, 4, 5, 10};
    } else {
        skinSizes = {-1};
    }

    const auto numbersSize = numbers.size();
    const auto skinSizesSize = skinSizes.size();

    // results are performance times and the (N,L,S,R) observables
    readdy::scalar timeForces[numbersSize][skinSizesSize];
    readdy::scalar timeIntegrator[numbersSize][skinSizesSize];
    readdy::scalar timeNeighborList[numbersSize][skinSizesSize];
    readdy::scalar timeReactions[numbersSize][skinSizesSize];
    readdy::scalar particleNumber[numbersSize][skinSizesSize];
    readdy::scalar systemSize[numbersSize][skinSizesSize];
    readdy::scalar relativeDisplacement[numbersSize][skinSizesSize];
    readdy::scalar reactivity[numbersSize][skinSizesSize];

    readdy::log::console()->set_level(spdlog::level::warn);
    for (auto i = 0; i < numbersSize; ++i) {
        for (auto j = 0; j < skinSizesSize; ++j) {
            readdy::log::error("at numbers {} / {} = {}, skins {} / {} = {}", i, numbersSize, numbers[i], j,
                               skinSizesSize, skinSizes[j]);
            std::map<std::string, readdy::scalar> factors;
            factors.emplace(std::make_pair(NUMBERS_FACTOR, numbers[i]));
            factors.emplace(std::make_pair(BOXLENGTHS_FACTOR, std::cbrt(numbers[i])));
            factors.emplace(std::make_pair(SKIN_FACTOR, skinSizes[j]));
            Scenario_t scenario(kernelName, factors);
            readdy::log::warn("this places us at {} particles", scenario.getParticleNumber());
            scenario.template perform<ReactionScheduler>(10, true);
            timeForces[i][j] = scenario.getTimeForces();
            timeIntegrator[i][j] = scenario.getTimeIntegrator();
            timeNeighborList[i][j] = scenario.getTimeNeighborlist();
            timeReactions[i][j] = scenario.getTimeReactions();
            particleNumber[i][j] = scenario.getParticleNumber();
            systemSize[i][j] = scenario.getSystemSize();
            relativeDisplacement[i][j] = scenario.getRelativeDisplacement();
            reactivity[i][j] = scenario.getReactivity();
        }
    }
    readdy::log::console()->set_level(spdlog::level::debug);

    {
        /** Write the result vectors to one dataset each, all in a single file */
        const std::string filename = "numbersSkinsConstDensity_" + Scenario_t::name + "_" + kernelName + ".h5";
        file_t file(filename, file_t::Action::CREATE, file_t::Flag::OVERWRITE);

        auto inputGroup = file.createGroup("/input");
        inputGroup.write(NUMBERS_FACTOR, numbers);
        inputGroup.write(SKIN_FACTOR, skinSizes);

        auto outputGroup = file.createGroup("/output");
        outputGroup.write("time_forces", {numbersSize, skinSizesSize}, &timeForces[0][0]);
        outputGroup.write("time_integrator", {numbersSize, skinSizesSize}, &timeIntegrator[0][0]);
        outputGroup.write("time_neighborlist", {numbersSize, skinSizesSize}, &timeNeighborList[0][0]);
        outputGroup.write("time_reactions", {numbersSize, skinSizesSize}, &timeReactions[0][0]);
        outputGroup.write("particle_number", {numbersSize, skinSizesSize}, &particleNumber[0][0]);
        outputGroup.write("system_size", {numbersSize, skinSizesSize}, &systemSize[0][0]);
        outputGroup.write("relative_displacement", {numbersSize, skinSizesSize}, &relativeDisplacement[0][0]);
        outputGroup.write("reactivity", {numbersSize, skinSizesSize}, &reactivity[0][0]);
    }
}

/*template<typename Scenario_t, typename ReactionScheduler=readdy::model::actions::reactions::Gillespie>
void scaleNumbersAndSkin_tmp(const std::string kernelName, bool reducedNumbers) {
    readdy::log::debug("using reduced numbers: {}", reducedNumbers);
    std::vector<readdy::scalar> numbers = {80};
    std::vector<readdy::scalar> skinSizes = {1.0};

    const auto numbersSize = numbers.size();
    const auto skinSizesSize = skinSizes.size();

    // results are performance times and the (N,L,S,R) observables
    readdy::scalar timeForces[numbersSize][skinSizesSize];
    readdy::scalar timeIntegrator[numbersSize][skinSizesSize];
    readdy::scalar timeNeighborList[numbersSize][skinSizesSize];
    readdy::scalar timeReactions[numbersSize][skinSizesSize];
    readdy::scalar particleNumber[numbersSize][skinSizesSize];
    readdy::scalar systemSize[numbersSize][skinSizesSize];
    readdy::scalar relativeDisplacement[numbersSize][skinSizesSize];
    readdy::scalar reactivity[numbersSize][skinSizesSize];

    for (auto i = 0; i < numbersSize; ++i) {
        for (auto j = 0; j < skinSizesSize; ++j) {
            readdy::log::error("at numbers {} / {} = {}, skins {} / {} = {}", i, numbersSize, numbers[i], j, skinSizesSize, skinSizes[j]);
            std::map<std::string, readdy::scalar> factors;
            factors.emplace(std::make_pair(NUMBERS_FACTOR, numbers[i]));
            factors.emplace(std::make_pair(BOXLENGTHS_FACTOR, std::cbrt(numbers[i])));
            factors.emplace(std::make_pair(SKIN_FACTOR, skinSizes[j]));
            Scenario_t scenario(kernelName, factors);
            readdy::log::warn("this places us at {} particles", scenario.getParticleNumber());
            scenario.template perform<ReactionScheduler>(20, true);
            timeForces[i][j] = scenario.getTimeForces();
            timeIntegrator[i][j] = scenario.getTimeIntegrator();
            timeNeighborList[i][j] = scenario.getTimeNeighborlist();
            timeReactions[i][j] = scenario.getTimeReactions();
            particleNumber[i][j] = scenario.getParticleNumber();
            systemSize[i][j] = scenario.getSystemSize();
            relativeDisplacement[i][j] = scenario.getRelativeDisplacement();
            reactivity[i][j] = scenario.getReactivity();
        }
    }

    {
        const std::string filename = "cpus=24.h5";
        file_t file(filename, file_t::Action::CREATE, file_t::Flag::OVERWRITE);

        auto inputGroup = file.createGroup("/input");
        inputGroup.write(NUMBERS_FACTOR, numbers);
        inputGroup.write(SKIN_FACTOR, skinSizes);

        auto outputGroup = file.createGroup("/output");
        outputGroup.write("time_forces", {numbersSize, skinSizesSize}, &timeForces[0][0]);
        outputGroup.write("time_integrator", {numbersSize, skinSizesSize}, &timeIntegrator[0][0]);
        outputGroup.write("time_neighborlist", {numbersSize, skinSizesSize}, &timeNeighborList[0][0]);
        outputGroup.write("time_reactions", {numbersSize, skinSizesSize}, &timeReactions[0][0]);
        outputGroup.write("particle_number", {numbersSize, skinSizesSize}, &particleNumber[0][0]);
        outputGroup.write("system_size", {numbersSize, skinSizesSize}, &systemSize[0][0]);
        outputGroup.write("relative_displacement", {numbersSize, skinSizesSize}, &relativeDisplacement[0][0]);
        outputGroup.write("reactivity", {numbersSize, skinSizesSize}, &reactivity[0][0]);
    }
}*/

TEST(TestPerformance, ReactiveCPU) {
    //scaleNumbersAndSkin<CollisiveUniformHomogeneous, readdy::model::actions::reactions::Gillespie>("SingleCPU", true);
    //scaleNumbersAndSkin<CollisiveUniformHomogeneous, readdy::model::actions::reactions::GillespieParallel>("CPU_Dense", false);
    //scaleNumbersAndSkin<ReactiveCollisiveUniformHomogeneous, readdy::model::actions::reactions::GillespieParallel>(
    //"CPU", false);
    readdy::Simulation sim;
    sim.setKernel("CPU");
    sim.setBoxSize(10, 10, 10);
    sim.setPeriodicBoundary({true, true, true});
    sim.registerParticleType("A", 1.0, 1.0);
    sim.registerHarmonicRepulsionPotential("A", "A", 1.0);
    for (int i = 0; i < 20; ++i) {
        sim.addParticle("A", 0, 0, 0);
    }
    readdy::log::console()->set_level(spdlog::level::debug);
    {
        using timer = readdy::util::RAIITimer;
        timer c(true, "timer");
        sim.runScheme(true).includeForces(false).withReactionScheduler("Gillespie")
                .evaluateObservables(false).configureAndRun(10, 1);
        //readdy::log::debug("Time {}", readdy::util::Timer::times().at("timer"));
    }

}

}
