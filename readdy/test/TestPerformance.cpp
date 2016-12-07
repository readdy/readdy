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

#include <hdf5.h>
#include <hdf5_hl.h>

#include <gtest/gtest.h>

#include <readdy/testing/Utils.h>
#include <readdy/common/Timer.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/Utils.h>

namespace {

/**
 * A performance scenario is an instance of a simulation that is run once for a few timesteps to record
 * the duration of force calculation, diffusion integration, neighborlist updating and performing reactions.
 */
class PerformanceScenario {
public:
    PerformanceScenario(const std::string kernelName) {
        kernel = readdy::plugin::KernelProvider::getInstance().create(kernelName);
        neighborList = kernel->createProgram<readdy::model::programs::UpdateNeighborList>();
    }

    template<typename ReactionScheduler=readdy::model::programs::reactions::Gillespie>
    void perform(const readdy::model::observables::time_step_type steps = 5) {
        this->configure();
        const auto result = runPerformanceTest<ReactionScheduler>(steps, false);
        timeForces = std::get<0>(result);
        timeIntegrator = std::get<1>(result);
        timeNeighborlist = std::get<2>(result);
        timeReactions = std::get<3>(result);
        hasPerformed = true;
    };

    double getTimeForces() const {
        if (hasPerformed) return timeForces;
        else throw std::runtime_error("scenario has not performed yet");
    }

    double getTimeIntegrator() const {
        if (hasPerformed) return timeIntegrator;
        else throw std::runtime_error("scenario has not performed yet");
    }

    double getTimeNeighborlist() const {
        if (hasPerformed) return timeNeighborlist;
        else throw std::runtime_error("scenario has not performed yet");
    }

    double getTimeReactions() const {
        if (hasPerformed) return timeReactions;
        else throw std::runtime_error("scenario has not performed yet");
    }

    double getParticleNumber() const { return particleNumber; }

    double getSystemSize() const { return systemSize; }

    double getRelativeDisplacement() const { return relativeDisplacement; }

    double getReactivity() const { return reactivity; }

protected:
    /**
     * Run a simulation for a given number of timesteps. All configuration is included in the kernel.
     * Which reaction scheduler is used is determined via template parameters
     *
     * @param steps number of timesteps that will be performed
     * @param verbose decides if to print information on context and performance
     * @return computation time per timestep for {forces, integrator, neighborlist, reactions}
     */
    template<typename ReactionScheduler=readdy::model::programs::reactions::Gillespie>
    std::array<double, 4>
    runPerformanceTest(readdy::model::observables::time_step_type steps, const bool verbose = false) {
        auto &&integrator = kernel->createProgram<readdy::model::programs::EulerBDIntegrator>();
        auto &&forces = kernel->createProgram<readdy::model::programs::CalculateForces>();
        auto &&reactionsProgram = kernel->createProgram<ReactionScheduler>();
        kernel->getKernelContext().configure(verbose);

        using timer = readdy::util::Timer;
        double timeForces = 0, timeIntegrator = 0, timeNeighborList = 0, timeReactions = 0;

        neighborList->execute();
        forces->execute();
        for (readdy::model::observables::time_step_type t = 0; t < steps; ++t) {
            if (verbose) {
                readdy::log::console()->debug("----------");
                readdy::log::console()->debug("t = {}", t);
            }
            {
                timer c("integrator", verbose);
                integrator->execute();
                timeIntegrator += c.getSeconds();
            }
            {
                timer c("neighbor list 1", verbose);
                neighborList->execute();
                timeNeighborList += c.getSeconds();
            }
            {
                timer c("reactions", verbose);
                reactionsProgram->execute();
                timeReactions += c.getSeconds();
            }
            {
                timer c("neighbor list 2", verbose);
                neighborList->execute();
                timeNeighborList += c.getSeconds();
            }
            {
                timer c("forces", verbose);
                forces->execute();
                timeForces += c.getSeconds();
            }
        }
        neighborList->setAction(readdy::model::programs::UpdateNeighborList::Action::clear);
        neighborList->execute();
        timeForces /= steps;
        timeIntegrator /= steps;
        timeNeighborList /= steps;
        timeReactions /= steps;
        if (verbose) {
            std::cout << "--------------------------------------------------------------" << std::endl;
            std::cout << "Average time for calculating forces: " << timeForces << std::endl;
            std::cout << "Average time for the integrator:     " << timeIntegrator << std::endl;
            std::cout << "Average time for the neighbor list:  " << timeNeighborList << std::endl;
            std::cout << "Average time for handling reactions: " << timeReactions << std::endl;
            std::cout << "--------------------------------------------------------------" << std::endl;
        }
        const std::array<double, 4> result = {timeForces, timeIntegrator, timeNeighborList, timeReactions};
        return result;
    }

    /** This is where the actual system gets set up: define reactions/potentials, add particles ... */
    virtual void configure() = 0;

    double timeForces, timeIntegrator, timeNeighborlist, timeReactions; // main result
    double particleNumber, systemSize, relativeDisplacement, reactivity;

protected:
    // the (N,L,S,R) observables defined above
    std::unique_ptr<readdy::model::Kernel> kernel;
    std::unique_ptr<readdy::model::programs::UpdateNeighborList> neighborList;
    bool hasPerformed = false;
};

/**
 * The following scenario shall simulate a purely reactive system. Particles are uniformly distributed and
 * particle-radii/reaction-radii are homogeneous, i.e. there is only one lengthscale involved, the reaction radius = 4.5.
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

    ReactiveUniformHomogeneous(const std::string kernelName, const std::map<std::string, double> factors)
            : PerformanceScenario(kernelName) {
        numberA = static_cast<unsigned long>(numberA * factors.at("numbers"));
        numberC = static_cast<unsigned long>(numberC * factors.at("numbers"));
        boxlength *= factors.at("boxlengths");
        neighborList->setSkinSize(factors.at("skin"));
        // set (N,L,S,R) observables
        particleNumber = 2 * numberA + numberC;
        systemSize = std::pow(boxlength, 3.) / std::pow(4.5, 3.);
        relativeDisplacement = std::sqrt(2 * diffusionA * timestep) / 4.5;
        reactivity = rateOn * timestep;
    };

    virtual void configure() override {
        auto &ctx = kernel->getKernelContext();
        ctx.setKBT(1.);
        ctx.setBoxSize(boxlength, boxlength, boxlength);
        ctx.setPeriodicBoundary(true, true, true);
        ctx.setDiffusionConstant("A", diffusionA);
        ctx.setDiffusionConstant("B", diffusionB);
        ctx.setDiffusionConstant("C", diffusionC);
        ctx.setTimeStep(timestep);

        kernel->registerReaction<readdy::model::reactions::Fusion>("A+B->C", "A", "B", "C", rateOn, 4.5);
        kernel->registerReaction<readdy::model::reactions::Fission>("C->A+B", "C", "A", "B", rateOff, 4.5);

        /** Distribute particles uniformly in box */
        auto &state = kernel->getKernelStateModel();
        std::vector<readdy::model::Particle> particles;
        particles.reserve(2*numberA + numberC);
        auto uniform = [this]() {
            return readdy::model::rnd::uniform_real<double, std::mt19937>(-0.5 * boxlength, 0.5 * boxlength);
        };
        const auto &typeMapping = ctx.getTypeMapping();
        for (auto i = 0; i < numberA; ++i) {
            readdy::model::Particle particleA(uniform(), uniform(), uniform(), typeMapping.at("A"));
            particles.push_back(particleA);
            readdy::model::Particle particleB(uniform(), uniform(), uniform(), typeMapping.at("B"));
            particles.push_back(particleB);
        }
        for (auto i = 0; i < numberC; ++i) {
            readdy::model::Particle particleC(uniform(), uniform(), uniform(), typeMapping.at("C"));
            particles.push_back(particleC);
        }
        state.addParticles(particles);
    }

private:
    unsigned long numberA = 500, numberC = 1800;
    double boxlength = 100., rateOn = 1e-3, rateOff = 5e-5, diffusionA = 0.14, diffusionB = 0.07, diffusionC = 0.06, timestep = 0.1;
};

const std::string ReactiveUniformHomogeneous::name = "ReactiveUniformHomogeneous";

template<typename scenario_t, typename ReactionScheduler=readdy::model::programs::reactions::Gillespie>
void scaleNumbersAndBoxsize(const std::string kernelName) {
    /** Base values will be multiplied by factors. numbers[i] and boxlength[i] factors for same i will conserve particle density */
    const std::vector<double> numbers = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 4., 5., 6., 7.};
    std::vector<double> boxlengths(numbers.size());
    {
        auto x = 0;
        std::generate(boxlengths.begin(), boxlengths.end(), [&numbers, &x] { return std::cbrt(numbers[x++]); });
    }
    bool supportsSkin = false;
    {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create(kernelName);
        auto program = kernel.get()->createProgram<readdy::model::programs::UpdateNeighborList>();
        supportsSkin = program->supportsSkin();
    }
    const auto numbersSize = numbers.size();
    const auto boxlengthsSize = boxlengths.size();

    // results are performance times and the (N,L,S,R) observables
    double timeForces[numbersSize][boxlengthsSize];
    double timeIntegrator[numbersSize][boxlengthsSize];
    double timeNeighborlist[numbersSize][boxlengthsSize];
    double timeReactions[numbersSize][boxlengthsSize];
    double particleNumber[numbersSize][boxlengthsSize];
    double systemSize[numbersSize][boxlengthsSize];
    double relativeDisplacement[numbersSize][boxlengthsSize];
    double reactivity[numbersSize][boxlengthsSize];
    double skinSize[numbersSize][boxlengthsSize];

    for (auto n = 0; n < numbers.size(); ++n) {
        for (auto l = 0; l < boxlengths.size(); ++l) {
            std::map<std::string, double> factors;
            factors.emplace(std::make_pair("numbers", numbers[n]));
            factors.emplace(std::make_pair("boxlengths", boxlengths[l]));
            scenario_t scenario(kernelName, factors);
            scenario.template perform<ReactionScheduler>(5);
            timeForces[n][l] = scenario.getTimeForces();
            timeIntegrator[n][l] = scenario.getTimeIntegrator();
            timeNeighborlist[n][l] = scenario.getTimeNeighborlist();
            timeReactions[n][l] = scenario.getTimeReactions();
            particleNumber[n][l] = scenario.getParticleNumber();
            systemSize[n][l] = scenario.getSystemSize();
            relativeDisplacement[n][l] = scenario.getRelativeDisplacement();
            reactivity[n][l] = scenario.getReactivity();
        }
    }

    /** Write the result vectors to one dataset each, all in a single file */
    const std::string filename = "numbersAndBoxlengths_" + scenario_t::name + "_" + kernelName + ".h5";
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t input = H5Gcreate(file, "/input", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dimsNumbers[1] = {numbersSize};
    hsize_t dimsBoxlengths[1] = {boxlengthsSize};
    H5LTmake_dataset_double(input, "numbers", 1, dimsNumbers, numbers.data());
    H5LTmake_dataset_double(input, "boxlengths", 1, dimsBoxlengths, boxlengths.data());

    hid_t output = H5Gcreate(file, "/output", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims_output[2] = {numbersSize, boxlengthsSize};
    H5LTmake_dataset_double(output, "time_forces", 2, dims_output, &timeForces[0][0]);
    H5LTmake_dataset_double(output, "time_integrator", 2, dims_output, &timeIntegrator[0][0]);
    H5LTmake_dataset_double(output, "time_neighborlist", 2, dims_output, &timeNeighborlist[0][0]);
    H5LTmake_dataset_double(output, "time_reactions", 2, dims_output, &timeReactions[0][0]);
    H5LTmake_dataset_double(output, "particle_number", 2, dims_output, &particleNumber[0][0]);
    H5LTmake_dataset_double(output, "system_size", 2, dims_output, &systemSize[0][0]);
    H5LTmake_dataset_double(output, "relative_displacement", 2, dims_output, &relativeDisplacement[0][0]);
    H5LTmake_dataset_double(output, "reactivity", 2, dims_output, &reactivity[0][0]);
    H5Fclose(file);
}

template<typename Scenario_t, typename ReactionScheduler=readdy::model::programs::reactions::Gillespie>
void scaleNumbersConstDensity(const std::string kernelName) {
    /** Base values will be multiplied by factors. numbers[i] and boxlength[i] factors for same i will conserve particle density */
    const std::vector<double> numbers = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 4., 5., 6., 7., 8., 9.,
                                         10., 12., 15., 20.};
    const auto numbersSize = numbers.size();

    // results are performance times and the (N,L,S,R) observables
    double timeForces[numbersSize];
    double timeIntegrator[numbersSize];
    double timeNeighborlist[numbersSize];
    double timeReactions[numbersSize];
    double particleNumber[numbersSize];
    double systemSize[numbersSize];
    double relativeDisplacement[numbersSize];
    double reactivity[numbersSize];

    for (auto n = 0; n < numbers.size(); ++n) {
        std::map<std::string, double> factors;
        factors.emplace(std::make_pair("numbers", numbers[n]));
        factors.emplace(std::make_pair("boxlengths", std::cbrt(numbers[n])));
        Scenario_t scenario(kernelName, factors);
        scenario.template perform<ReactionScheduler>(20);
        timeForces[n] = scenario.getTimeForces();
        timeIntegrator[n] = scenario.getTimeIntegrator();
        timeNeighborlist[n] = scenario.getTimeNeighborlist();
        timeReactions[n] = scenario.getTimeReactions();
        particleNumber[n] = scenario.getParticleNumber();
        systemSize[n] = scenario.getSystemSize();
        relativeDisplacement[n] = scenario.getRelativeDisplacement();
        reactivity[n] = scenario.getReactivity();
    }

    /** Write the result vectors to one dataset each, all in a single file */
    const std::string filename = "numbersConstDensity_" + Scenario_t::name + "_" + kernelName + ".h5";
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t input = H5Gcreate(file, "/input", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dimsNumbers[1] = {numbersSize};
    H5LTmake_dataset_double(input, "numbers", 1, dimsNumbers, numbers.data());

    hid_t output = H5Gcreate(file, "/output", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims_output[1] = {numbersSize};
    H5LTmake_dataset_double(output, "time_forces", 1, dims_output, &timeForces[0]);
    H5LTmake_dataset_double(output, "time_integrator", 1, dims_output, &timeIntegrator[0]);
    H5LTmake_dataset_double(output, "time_neighborlist", 1, dims_output, &timeNeighborlist[0]);
    H5LTmake_dataset_double(output, "time_reactions", 1, dims_output, &timeReactions[0]);
    H5LTmake_dataset_double(output, "particle_number", 1, dims_output, &particleNumber[0]);
    H5LTmake_dataset_double(output, "system_size", 1, dims_output, &systemSize[0]);
    H5LTmake_dataset_double(output, "relative_displacement", 1, dims_output, &relativeDisplacement[0]);
    H5LTmake_dataset_double(output, "reactivity", 1, dims_output, &reactivity[0]);
    H5Fclose(file);
}

TEST(TestPerformance, ReactiveCPU) {
    scaleNumbersConstDensity<ReactiveUniformHomogeneous, readdy::model::programs::reactions::Gillespie>("CPU");
}

}