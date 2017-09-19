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


#include <spdlog/spdlog.h>
#include <readdy/api/Simulation.h>
#include <readdy/common/range.h>
#include <readdy/common/Timer.h>

/**
 *
 *
 * @file DiffusionBenchmark.cpp
 * @brief
 * @author clonker
 * @author chrisfroe
 * @date 7/20/17
 */

int main(int argc, char **argv) {
    auto console = spdlog::stdout_color_mt("console");
    console->set_level(spdlog::level::debug);
    console->set_pattern("[          ] [%Y-%m-%d %H:%M:%S] [%t] [%l] %v");
    readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory("readdy/readdy_plugins");

    std::size_t n_timesteps = 1000;
    readdy::scalar timestep = .1;

    readdy::scalar box_length = 450.;
    int n_particles = 1000;
    readdy::scalar diffusion_coefficient = 1.;
    readdy::scalar particle_radius = 1.25;
    readdy::scalar force_constant = 100.;

    readdy::Simulation simulation;
    simulation.setKernel("CPU");
    simulation.setPeriodicBoundary({true, true, true});
    simulation.setBoxSize(readdy::Vec3(box_length, box_length, box_length));
    simulation.registerParticleType("A", diffusion_coefficient, particle_radius);
    simulation.registerHarmonicRepulsionPotential("A", "A", force_constant);

    for(auto _ : readdy::util::range<int>(0, n_particles)) {
        simulation.addParticle("A",
                               box_length * readdy::model::rnd::uniform_real(0., 1.) - .5 * box_length,
                               box_length * readdy::model::rnd::uniform_real(0., 1.) - .5 * box_length,
                               box_length * readdy::model::rnd::uniform_real(0., 1.) - .5 * box_length);
    }

    using timer = readdy::util::Timer;

    {
        //timer clock {"timer"};
        simulation.runScheme<readdy::api::ReaDDyScheme>(true).withSkinSize(50.).configureAndRun(n_timesteps, timestep);
    }

    return 0;
}