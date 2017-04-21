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
 * << detailed description >>
 *
 * @file TestBug.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <gtest/gtest.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/api/SimulationScheme.h>

namespace  {

TEST(TestBug, DiffusionAndReaction) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
    kernel->getKernelContext().particle_types().add("A", 0.05, 1.0);
    kernel->getKernelContext().particle_types().add("F", 0.0, 1.0);
    kernel->getKernelContext().particle_types().add("V", 0.0, 1.0);
    kernel->getKernelContext().setPeriodicBoundary(true, true, true);
    kernel->getKernelContext().setBoxSize(100, 10, 10);

    const double weightF = 0.;
    const double weightA = 1.;
    kernel->registerReaction<readdy::model::reactions::Fusion>("F+A->F", "A", "F", "F", .1, 2.0, weightF, weightA);
    //kernel->registerReaction<readdy::model::reactions::Fusion>("F+A->F2", "V", "V", "V", .1, 2.0, weightF, weightA);

    auto n3 = readdy::model::rnd::normal3<>;
    // 120 F particles
    for (std::size_t i = 0; i < 100; ++i) {
        kernel->addParticle("F", n3(0., 1.));
        kernel->addParticle("A", n3(0., 1.));
    }

    auto obs = kernel->createObservable<readdy::model::observables::NParticles>(1, std::vector<std::string>({"F", "A"}));
    obs->setCallback(
            [&](const readdy::model::observables::NParticles::result_t &result) {
                log::warn("#F={}, #A={}", result[0], result[1]);
            }
    );
    auto connection = kernel->connectObservable(obs.get());

    {
        auto conf = readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>(kernel.get(), true);
        conf.withReactionScheduler<readdy::model::actions::reactions::Gillespie>();
        conf.withSkinSize(.1);
        conf.configureAndRun(.01, 10);
    }
}

TEST(TestBug, Diffusion) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
    auto& context = kernel->getKernelContext();
    context.particle_types().add("A", 0.05, 1.0);
    context.particle_types().add("F", 0.0, 1.0);
    context.particle_types().add("V", 0.0, 1.0);
    // just to have a cutoff
    kernel->registerReaction<readdy::model::reactions::Fusion>("test", "V", "V", "V", .1, 2.0);
    context.setPeriodicBoundary(true, true, true);
    context.setBoxSize(100, 10, 10);

    auto n3 = readdy::model::rnd::normal3<>;
    // 120 F particles
    for (std::size_t i = 0; i < 100; ++i) {
        kernel->addParticle("F", n3(0., 1.));
        kernel->addParticle("A", n3(0., 1.));
    }
    auto obs = kernel->createObservable<readdy::model::observables::NParticles>(1);
    obs->setCallback(
            [&](const readdy::model::observables::NParticles::result_t &result) {
                bool wrong_i_j = false, wrong_j_i = false;
                auto d2 = context.getDistSquaredFun();
                const auto neighbor_list = kernel->getCPUKernelStateModel().getNeighborList();
                const auto& data = *kernel->getCPUKernelStateModel().getParticleData();
                std::size_t i = 0;
                for(auto it_i = data.begin(); it_i != data.end(); ++it_i, ++i) {
                    std::size_t j = 0;
                    for(auto it_j = data.begin(); it_j != data.end(); ++it_j, ++j) {
                        if(it_i != it_j && std::sqrt(d2(it_i->position(), it_j->position())) < 2.0) {
                            auto neigh_i = neighbor_list->find_neighbors(i);
                            auto neigh_j = neighbor_list->find_neighbors(j);
                            if(std::find(neigh_i.begin(), neigh_i.end(), j) == neigh_i.end()) {
                                wrong_i_j = true;
                            }
                            if(std::find(neigh_j.begin(), neigh_j.end(), i) == neigh_j.end()) {
                                wrong_j_i = true;
                            }
                        }
                    }
                }
                if(wrong_i_j || wrong_j_i) {
                    log::critical("WRONG: i->j {}, j->i {}", wrong_i_j, wrong_j_i);
                } else {
                    log::error("this worked");
                }
            }
    );
    auto connection = kernel->connectObservable(obs.get());
    {
        auto conf = readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>(kernel.get(), true);
        conf.withReactionScheduler<readdy::model::actions::reactions::Gillespie>();
        conf.withSkinSize(.1);
        conf.configureAndRun(.01, 4);
    }
}

TEST(TestBug, MarkDirty) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
    auto& context = kernel->getKernelContext();
    context.particle_types().add("A", 0.05, 1.0);
    context.particle_types().add("F", 0.0, 1.0);
    context.particle_types().add("V", 0.0, 1.0);
    // just to have a cutoff
    kernel->registerReaction<readdy::model::reactions::Fusion>("test", "V", "V", "V", .1, 2.0);
    context.setPeriodicBoundary(true, true, true);
    context.setBoxSize(100, 10, 10);
    kernel->addParticle("A", model::Vec3());

    {
        auto& stateModel = kernel->getCPUKernelStateModel();
        auto neighborList = stateModel.getNeighborList();
        neighborList->setSkinSize(1.0);
        neighborList->create();
        auto& data = *stateModel.getParticleData();
        neighborList->displace(data.begin(), {1.0, 0, 0});
    }
}

}