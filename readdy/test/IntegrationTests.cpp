//
// Created by mho on 6/17/19.
//

#include <catch2/catch.hpp>

#include <readdy/readdy.h>
#include <readdy/common/FloatingPoints.h>

#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

using namespace readdytesting::kernel;

static constexpr const char *const REACTION_HANDLERS[] = {
        "UncontrolledApproximation", "Gillespie", "DetailedBalance"
};

TEMPLATE_TEST_CASE("Reaction handlers integration.", "[!hide][integration]", SingleCPU, CPU) {
    auto kernel = create<TestType>();
    auto &ctx = kernel->context();
    for (const auto &handler : REACTION_HANDLERS) {
        SECTION(handler) {
            if (kernel->name() == "CPU" && std::string(handler) == "DetailedBalance") {
                INFO("Skipping tests for the combination (CPU, DetailedBalance)");
                continue;
            }

            SECTION("Michaelis Menten") {
                /**
                * Since comparing the value of a stochastic process (number of particles over time) is not well
                * suited for testing, here integrate the stochastic process over the simulated time,
                * giving a quite robust observable. This can be compared to an analytic value from ODE kinetics.
                */
                namespace rnd = readdy::model::rnd;

                readdy::scalar length {10.};
                ctx.boxSize() = {{length, length, length}};
                ctx.periodicBoundaryConditions() = {true, true, true};

                ctx.particleTypes().add("E", 0.5);
                ctx.particleTypes().add("S", 0.5);
                ctx.particleTypes().add("ES", 0.5);
                ctx.particleTypes().add("P", 0.5);

                ctx.reactions().add("fwd: E +(1.0) S -> ES", 0.0023417691399750494);
                ctx.reactions().add("back: ES -> E +(1.0) S", 1.);
                ctx.reactions().add("prod: ES -> E +(1.0) P", 1.);

                readdy::scalar dt {8e-3};
                auto &&reactions = kernel->actions().createReactionScheduler(handler, dt);
                auto &&bd = kernel->actions().eulerBDIntegrator(dt);
                auto &&initNeighborList = kernel->actions().createNeighborList(1.0);
                auto &&neighborList = kernel->actions().updateNeighborList();

                std::size_t nParticlesE {70};
                std::size_t nParticlesS {700};

                for (int i = 0; i < nParticlesE; ++i) {
                    readdy::Vec3 pos{rnd::uniform_real() * length - 0.5 * length,
                                     rnd::uniform_real() * length - 0.5 * length,
                                     rnd::uniform_real() * length - 0.5 * length};
                    kernel->addParticle("E", pos);
                }

                for (int i = 0; i < nParticlesS; ++i) {
                    readdy::Vec3 pos{rnd::uniform_real() * length - 0.5 * length,
                                     rnd::uniform_real() * length - 0.5 * length,
                                     rnd::uniform_real() * length - 0.5 * length};
                    kernel->addParticle("S", pos);
                }

                std::size_t nSteps {37500};

                std::vector<std::string> typesToCount {{"S"}};
                auto &&obs = kernel->observe().nParticles(1, typesToCount);
                readdy::scalar meanS {0.};

                obs->callback() = [&](const readdy::model::observables::NParticles::result_type &result) {
                    meanS += static_cast<readdy::scalar>(result[0]);
                };

                auto &&connection = kernel->connectObservable(obs.get());

                initNeighborList->perform();
                neighborList->perform();
                kernel->evaluateObservables(0);
                auto t1 = std::chrono::high_resolution_clock::now();
                for (readdy::TimeStep t = 1; t < nSteps+1; t++) {
                    if (t == (nSteps+1)/100) {
                        auto t2 = std::chrono::high_resolution_clock::now();
                        auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                        readdy::log::warn("1% done in {} ms, ETA is {} seconds", dur, dur / 10.);
                    }
                    bd->perform();
                    neighborList->perform();
                    reactions->perform();
                    neighborList->perform();
                    kernel->evaluateObservables(t);
                }

                meanS /= static_cast<readdy::scalar>(nSteps+1);

                CHECK(meanS < 665.2 + 0.03 * 665.2);
                CHECK(meanS > 665.2 - 0.03 * 665.2);
            }
        }
    }
}

TEMPLATE_TEST_CASE("Chain-decay integration test", "[!hide][integration]", SingleCPU, CPU) {
    auto kernel = create<TestType>();
    auto &ctx = kernel->context();
    ctx.boxSize() = {{150, 150, 150}};
    ctx.periodicBoundaryConditions() = {{true, true, true}};

    ctx.particleTypes().add("Decay", 1.);
    ctx.particleTypes().addTopologyType("T", 1.);
    ctx.particleTypes().add("whatever", 1.);
    ctx.reactions().add("decay: Decay ->", 1000.);
    ctx.reactions().add("whatever: whatever ->", 1000.);
    ctx.potentials().addBox("Decay", 10, {-70, -70, -70}, {130, 130, 130});
    ctx.potentials().addBox("T", 10, {-70, -70, -70}, {130, 130, 130});

    ctx.topologyRegistry().configureBondPotential("T", "T", {20, 2});
    ctx.topologyRegistry().addType("unstable");
    // decay reaction
    auto reactionFunction = [&](readdy::model::top::GraphTopology &top) {
        readdy::model::top::reactions::Recipe recipe (top);
        auto& vertices = top.graph().vertices();
        auto current_n_vertices = vertices.size();
        auto vIdx = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices - 1));
        auto v = vertices.begin();
        std::advance(v, vIdx);
        recipe.separateVertex(v);
        recipe.changeParticleType(v, "Decay");
        return recipe;
    };
    auto rateFunction = [](const readdy::model::top::GraphTopology &top) {
        return .1;
    };
    readdy::model::top::reactions::StructuralTopologyReaction reaction {"rr", reactionFunction, rateFunction};
    ctx.topologyRegistry().addStructuralReaction("unstable", reaction);

    std::vector<readdy::model::TopologyParticle> topologyParticles;
    {
        topologyParticles.reserve(70);
        for (std::size_t i = 0; i < 70; ++i) {
            const auto id = ctx.particleTypes().idOf("T");
            topologyParticles.emplace_back(-5 + i * 10. / static_cast<readdy::scalar>(70), 0, 0, id);
        }
    }

    for(auto i = 0U; i < 50; ++i) {
        kernel->stateModel().addParticle(readdy::model::Particle(0, 0, 0, ctx.particleTypes().idOf("whatever")));
    }

    auto topology = kernel->stateModel().addTopology(ctx.topologyRegistry().idOf("unstable"), topologyParticles);
    {
        auto it = topology->graph().vertices().begin();
        auto it2 = ++topology->graph().vertices().begin();
        while(it2 != topology->graph().vertices().end()) {
            topology->graph().addEdge(it, it2);
            std::advance(it, 1);
            std::advance(it2, 1);
        }
    }

    auto topsobs = kernel->observe().topologies(1);
    auto connection = kernel->connectObservable(topsobs.get());
    topsobs->callback() = [&](const readdy::model::observables::Topologies::result_type &result) {
        auto tops = kernel->stateModel().getTopologies();
        REQUIRE(result.size() == tops.size());
        for(std::size_t i = 0; i < tops.size(); ++i) {
            auto top = tops.at(i);
            auto record = result.at(i);

            auto edges = top->graph().edges();
            REQUIRE(edges.size() == record.edges.size());

            auto particles = top->fetchParticles();

            REQUIRE(record.particleIndices.size() == particles.size());
            for(std::size_t j = 0; j < record.particleIndices.size(); ++j) {
                auto topParticles = top->getParticles();
                kernel->stateModel().toDenseParticleIndices(topParticles.begin(), topParticles.end());
                REQUIRE(std::find(topParticles.begin(), topParticles.end(),
                                  record.particleIndices.at(j)) != topParticles.end());
                REQUIRE(particles.at(j).type() == ctx.particleTypes().idOf("T"));
            }

            for(const auto &edge : record.edges) {
                auto it = std::find_if(edges.begin(), edges.end(), [&](const auto &e) {
                    auto p1Idx = std::get<0>(e)->particleIndex;
                    auto p2Idx = std::get<1>(e)->particleIndex;
                    return (p1Idx == std::get<0>(edge) && p2Idx == std::get<1>(edge))
                           || (p1Idx == std::get<1>(edge) && p2Idx == std::get<0>(edge));
                });
                REQUIRE(it != edges.end());
            }

        }
    };

    {
        auto integrator = kernel->actions().createIntegrator("EulerBDIntegrator", 1e-2);
        auto forces = kernel->actions().calculateForces();
        auto topReactions = kernel->actions().evaluateTopologyReactions(1e-2);
        auto reactions = kernel->actions().uncontrolledApproximation(1e-2);

        std::size_t time = 0;
        std::size_t n_time_steps = 10000;

        kernel->initialize();

        forces->perform();
        kernel->evaluateObservables(time);
        for(time = 1; time < n_time_steps; ++time) {
            integrator->perform();
            topReactions->perform();
            forces->perform();
            reactions->perform();
            kernel->evaluateObservables(time);
            for(auto topPtr : kernel->stateModel().getTopologies()) {
                // check that all topologies are just containing T particles and their edges are also fine
                for(const auto &p : topPtr->fetchParticles()) {
                    REQUIRE(p.type() == ctx.particleTypes().idOf("T"));
                }
                for(auto edge : topPtr->graph().edges()) {
                    auto v1 = std::get<0>(edge);
                    auto v2 = std::get<1>(edge);
                    REQUIRE(topPtr->particleForVertex(v1).type() == ctx.particleTypes().idOf("T"));
                    REQUIRE(topPtr->particleForVertex(v2).type() == ctx.particleTypes().idOf("T"));
                }
            }
        }
    }
}

TEMPLATE_TEST_CASE("Attach particle to topology", "[!hide][integration]", SingleCPU, CPU) {
    readdy::model::Context ctx;

    ctx.periodicBoundaryConditions() = {{true, true, true}};
    ctx.topologyRegistry().addType("TA");
    ctx.boxSize() = {{15, 15, 15}};
    ctx.particleTypes().add("middle", 0., readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("end", 0., readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("A", 0.);
    ctx.topologyRegistry().configureBondPotential("middle", "middle", {.00000001, 1});
    ctx.topologyRegistry().configureBondPotential("middle", "end", {.00000001, 1});
    ctx.topologyRegistry().configureBondPotential("end", "end", {.00000001, 1});
    // register attach reaction that transforms (end, A) -> (middle, end)
    ctx.topologyRegistry().addSpatialReaction("attach: TA (end) + (A) -> TA (middle--end)",
                                                               1e10, 1.5);

    readdy::Simulation simulation(create<TestType>(), ctx);

    auto top = simulation.addTopology("TA", {simulation.createTopologyParticle("end", {-1., 0., 0.}),
                                             simulation.createTopologyParticle("middle", {0., 0., 0.}),
                                             simulation.createTopologyParticle("end", {1., 0., 0.})});
    {
        auto it = top->graph().vertices().begin();
        auto it2 = std::next(top->graph().vertices().begin());
        top->graph().addEdge(it, it2);
        ++it; ++it2;
        top->graph().addEdge(it, it2);
    }

    simulation.addParticle("A", -2., 0, 0);
    simulation.addParticle("A", -3., 0, 0);
    simulation.addParticle("A", -4., 0, 0);
    simulation.addParticle("A", 2., 0, 0);
    simulation.addParticle("A", 3., 0, 0);
    simulation.addParticle("A", 4., 0, 0);

    simulation.run(6, 1.);

    const auto& type_registry = simulation.context().particleTypes();

    REQUIRE(simulation.context().topologyRegistry().isSpatialReactionType("A"));
    REQUIRE(simulation.context().topologyRegistry().isSpatialReactionType("end"));
    REQUIRE_FALSE(simulation.context().topologyRegistry().isSpatialReactionType("middle"));
    REQUIRE(simulation.context().calculateMaxCutoff() == 1.5);

    REQUIRE(simulation.currentTopologies().size() == 1);
    auto chainTop = simulation.currentTopologies().at(0);
    REQUIRE(chainTop->getNParticles() == 3 /*original topology particles*/ + 6 /*attached particles*/);

    auto top_particles = simulation.stateModel().getParticlesForTopology(*chainTop);

    bool foundEndVertex {false};
    // check that graph is indeed linear
    for(std::size_t idx = 0; idx < chainTop->graph().vertices().size() && !foundEndVertex; ++idx) {
        auto prev_neighbor = std::next(chainTop->graph().vertices().begin(), idx);
        const auto& v_end = *prev_neighbor;
        if(v_end.particleType() == type_registry.idOf("end")) {
            foundEndVertex = true;

            REQUIRE(v_end.neighbors().size() == 1);
            REQUIRE(top_particles.at(idx).type() == type_registry.idOf("end"));

            using flouble = readdy::fp::FloatingPoint<readdy::scalar>;
            flouble x_end (top_particles.at(idx).pos().x);
            flouble y_end (top_particles.at(idx).pos().y);
            flouble z_end (top_particles.at(idx).pos().z);
            // the end particle of our topology sausage should be either at x=4 or x=-4
            REQUIRE((x_end.AlmostEquals(flouble{4.}) || x_end.AlmostEquals(flouble{-4.})));
            REQUIRE(y_end.AlmostEquals(flouble{0.})); // no diffusion going on
            REQUIRE(z_end.AlmostEquals(flouble{0.})); // no diffusion going on

            auto factor = x_end.AlmostEquals(flouble{4.}) ? 1. : -1.;

            // walk along topology sausage, check end particles are always at +-4, the other ones are of type middle
            auto next_neighbor = v_end.neighbors().at(0);
            std::size_t i = 0;
            while(next_neighbor->particleType() != type_registry.idOf("end") && i < 20 /*no endless loop*/) {
                auto next_idx = std::distance(chainTop->graph().vertices().begin(), next_neighbor);
                const auto& next_particle = top_particles.at(static_cast<std::size_t>(next_idx));
                auto predicted_pos = factor*4. - factor*(i+1)*1.;
                auto actual_pos = next_particle.pos().x;
                REQUIRE((flouble(actual_pos).AlmostEquals(flouble(predicted_pos))));
                REQUIRE((flouble(next_particle.pos().y)).AlmostEquals(flouble(0.)));
                REQUIRE((flouble(next_particle.pos().z)).AlmostEquals(flouble(0.)));
                if(next_neighbor->particleType() == type_registry.idOf("middle")) {
                    REQUIRE(next_neighbor->neighbors().size() == 2);
                    if(next_neighbor->neighbors().at(0) == prev_neighbor) {
                        prev_neighbor = next_neighbor;
                        next_neighbor = next_neighbor->neighbors().at(1);
                    } else {
                        prev_neighbor = next_neighbor;
                        next_neighbor = next_neighbor->neighbors().at(0);
                    }

                } else {
                    REQUIRE(next_neighbor->neighbors().size() == 1);
                }

                ++i;
            }

            REQUIRE(i == 3+6-1-1);
        }
    }
    REQUIRE(foundEndVertex);
}

TEMPLATE_TEST_CASE("Break bonds due to pulling", "[!hide][breakbonds][integration]", SingleCPU, CPU) {
    GIVEN("A linear polymer with breakable bonds") {
        // bond force constants are adjusted such that spontaneous breaking is unlikely,
        // i.e. RMSD of bond when fluctuating is smaller than required when breaking threshold
        // i.e. only the external pulling should break the bond
        // RMSD = 1/sqrt(2 beta force-const)
        // Stability threshold is 'dt D beta force-const < 1'
        // When starting in rod-like configuration, there is an entropic force
        // that collapses the polymer. To avoid that, put the polymer in a long thin tube like potential
        auto kernel = readdytesting::kernel::create<TestType>();
        auto &ctx = kernel->context();
        auto &stateModel = kernel->stateModel();
        auto &topReg = ctx.topologyRegistry();
        auto &types = ctx.particleTypes();

        ctx.kBT() = 0.01; // low temperature for sharp distribution
        ctx.periodicBoundaryConditions() = {{true, true, true}};
        topReg.addType("T1");
        topReg.addType("T2");
        ctx.boxSize() = {{20, 10, 10}};
        types.add("head", 0.1, readdy::model::particleflavor::TOPOLOGY);
        types.add("A", 0.1, readdy::model::particleflavor::TOPOLOGY);
        types.add("tail", 0.1, readdy::model::particleflavor::TOPOLOGY);
        topReg.configureBondPotential("head", "A", {10., 2});
        topReg.configureBondPotential("A", "A", {10., 4});
        topReg.configureBondPotential("A", "tail", {10., 2});
        ctx.potentials().addCylinder("A", 100., {0.,0.,0.}, {1., 0., 0.}, 0.01, true);
        ctx.potentials().addCylinder("head", 100., {0.,0.,0.}, {1., 0., 0.}, 0.01, true);
        ctx.potentials().addCylinder("tail", 100., {0.,0.,0.}, {1., 0., 0.}, 0.01, true);

        ctx.potentials().addBox("head", 10., {-4., -12.5, -12.5}, {0.000001, 25., 25.});

        std::vector<readdy::model::TopologyParticle> particles{
            {-4., 0., 0., types.idOf("head")},
            {-2., 0., 0., types.idOf("A")},
            {2., 0., 0., types.idOf("A")},
            {4., 0., 0., types.idOf("tail")},
        };
        REQUIRE(particles.size() == 4);

        auto graphTop = stateModel.addTopology(topReg.idOf("T1"), particles);
        {
            auto &graph = graphTop->graph();
            for (std::size_t i = 0; i < 3; ++i) {
                graph.addEdgeBetweenParticles(i, i + 1);
            }
        }

        readdy::model::actions::top::BreakConfig breakConfig;
        breakConfig.addBreakablePair(types.idOf("A"), types.idOf("A"), 1.0, 1.0);

        readdy::scalar timeStep = 0.0005;
        std::size_t nSteps = 100000;
        auto diffusion = kernel->actions().eulerBDIntegrator(timeStep);
        auto forces = kernel->actions().calculateForces();
        auto breakingBonds = kernel->actions().breakBonds(timeStep, breakConfig);

        WHEN("an external potential pulls the tail particle") {
            ctx.potentials().addBox("tail", 100., {+8, -12.5, -12.5}, {0.000001, 25., 25.});

            kernel->initialize();
            forces->perform();
            for (std::size_t t = 1; t < nSteps+1; t++) {
                diffusion->perform();
                breakingBonds->perform();
                forces->perform();
            }

            THEN("some bond will break because the energy threshold is exceeded") {
                auto topsAfter = stateModel.getTopologies();
                REQUIRE_FALSE(topsAfter.empty());
                REQUIRE(topsAfter.size() > 1);
                // all sub topologies are linear chains (each vertex has at most 2 neighbors)
                for (const auto &top : topsAfter) {
                    const auto &graph = top->graph();
                    const auto &vertices = graph.vertices();
                    for (const auto &v : vertices) {
                        CHECK(v.neighbors().size() <= 2);
                    }
                }
            }
        }

        WHEN("there is no potential pulling on the chain") {
            kernel->initialize();
            forces->perform();
            for (std::size_t t = 1; t < nSteps+1; t++) {
                diffusion->perform();
                breakingBonds->perform();
                forces->perform();
            }

            THEN("No bond will break spontaneously") {
                auto topsAfter = stateModel.getTopologies();
                REQUIRE_FALSE(topsAfter.empty());
                REQUIRE(topsAfter.size() == 1);
                const auto &graph = topsAfter.at(0)->graph();
                const auto &vertices = graph.vertices();
                for (const auto &v : vertices) {
                    CHECK(v.neighbors().size() <= 2);
                    CHECK(!v.neighbors().empty());
                }
            }
        }
    }
}
