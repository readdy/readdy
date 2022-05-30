//
// Created by mho on 6/17/19.
//

#include <catch2/catch_template_test_macros.hpp>

#include <readdy/readdy.h>
#include <readdy/common/FloatingPoints.h>

#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/common/boundary_condition_operations.h>

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

                readdy::scalar length{10.};
                ctx.boxSize() = {{length, length, length}};
                ctx.periodicBoundaryConditions() = {true, true, true};

                ctx.particleTypes().add("E", 0.5);
                ctx.particleTypes().add("S", 0.5);
                ctx.particleTypes().add("ES", 0.5);
                ctx.particleTypes().add("P", 0.5);

                ctx.reactions().add("fwd: E +(1.0) S -> ES", 0.0023417691399750494);
                ctx.reactions().add("back: ES -> E +(1.0) S", 1.);
                ctx.reactions().add("prod: ES -> E +(1.0) P", 1.);

                readdy::scalar dt{8e-3};
                auto &&reactions = kernel->actions().createReactionScheduler(handler, dt);
                auto &&bd = kernel->actions().eulerBDIntegrator(dt);
                auto &&initNeighborList = kernel->actions().createNeighborList(1.0);
                auto &&neighborList = kernel->actions().updateNeighborList();

                std::size_t nParticlesE{70};
                std::size_t nParticlesS{700};

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

                std::size_t nSteps{37500};

                std::vector<std::string> typesToCount{{"S"}};
                auto &&obs = kernel->observe().nParticles(1, typesToCount);
                readdy::scalar meanS{0.};

                obs->setCallback([&](const readdy::model::observables::NParticles::result_type &result) {
                    meanS += static_cast<readdy::scalar>(result[0]);
                });

                auto &&connection = kernel->connectObservable(obs.get());

                initNeighborList->perform();
                neighborList->perform();
                kernel->evaluateObservables(0);
                auto t1 = std::chrono::high_resolution_clock::now();
                for (readdy::TimeStep t = 1; t < nSteps + 1; t++) {
                    if (t == (nSteps + 1) / 100) {
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

                meanS /= static_cast<readdy::scalar>(nSteps + 1);

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
        readdy::model::top::reactions::Recipe recipe(top);
        auto &vertices = top.graph().vertices();
        auto current_n_vertices = vertices.size();
        auto vIdx = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices - 1));
        auto v = vertices.begin();
        std::advance(v, vIdx);
        recipe.separateVertex(v.persistent_index());
        recipe.changeParticleType(v.persistent_index(), "Decay");
        return recipe;
    };
    auto rateFunction = [](const readdy::model::top::GraphTopology &top) {
        return .1;
    };
    readdy::model::top::reactions::StructuralTopologyReaction reaction{"rr", reactionFunction, rateFunction};
    ctx.topologyRegistry().addStructuralReaction("unstable", reaction);

    std::vector<readdy::model::Particle> topologyParticles;
    {
        topologyParticles.reserve(70);
        for (std::size_t i = 0; i < 70; ++i) {
            const auto id = ctx.particleTypes().idOf("T");
            topologyParticles.emplace_back(-5 + i * 10. / static_cast<readdy::scalar>(70), 0, 0, id);
        }
    }

    for (auto i = 0U; i < 50; ++i) {
        kernel->stateModel().addParticle(readdy::model::Particle(0, 0, 0, ctx.particleTypes().idOf("whatever")));
    }

    auto topology = kernel->stateModel().addTopology(ctx.topologyRegistry().idOf("unstable"), topologyParticles);
    {
        auto it = topology->graph().vertices().begin();
        auto it2 = ++topology->graph().vertices().begin();
        while (it2 != topology->graph().vertices().end()) {
            topology->addEdge(it.persistent_index(), it2.persistent_index());
            std::advance(it, 1);
            std::advance(it2, 1);
        }
    }

    auto topsobs = kernel->observe().topologies(1);
    auto connection = kernel->connectObservable(topsobs.get());
    topsobs->setCallback([&](const readdy::model::observables::Topologies::result_type &result) {
        auto tops = kernel->stateModel().getTopologies();
        REQUIRE(result.size() == tops.size());
        for (std::size_t i = 0; i < tops.size(); ++i) {
            auto top = tops.at(i);
            auto record = result.at(i);

            auto edges = top->graph().edges();
            REQUIRE(edges.size() == record.edges.size());

            auto particles = top->fetchParticles();

            REQUIRE(record.particleIndices.size() == particles.size());
            auto topParticlesDense = top->particleIndices();
            kernel->stateModel().toDenseParticleIndices(topParticlesDense.begin(), topParticlesDense.end());
            for (std::size_t j = 0; j < record.particleIndices.size(); ++j) {
                REQUIRE(std::find(topParticlesDense.begin(), topParticlesDense.end(),
                                  record.particleIndices.at(j)) != topParticlesDense.end());
                REQUIRE(particles.at(j).type() == ctx.particleTypes().idOf("T"));
            }
            for (const auto &edge : record.edges) {
                // topology-local indices
                auto[recordV1, recordV2] = edge;
                // map to particle indices in record space
                auto p1Idx = record.particleIndices.at(recordV1);
                auto p2Idx = record.particleIndices.at(recordV2);
                // as we already tested that record and topology are consistent wrt particle indices
                // we can pick the "argwhere" of the dense topology indices to find the index in the
                // topParticleDense corresponding to our particle
                auto realP1IdxLocal = std::distance(
                        topParticlesDense.begin(), std::find(topParticlesDense.begin(),
                                                             topParticlesDense.end(), p1Idx));
                auto realP2IdxLocal = std::distance(topParticlesDense.begin(),
                                                    std::find(topParticlesDense.begin(),
                                                              topParticlesDense.end(),
                                                              p2Idx));
                // that index can be mapped back to space with blanks like this
                auto v1Ix = top->vertexIndexForParticle(top->particleIndices().at(realP1IdxLocal));
                auto v2Ix = top->vertexIndexForParticle(top->particleIndices().at(realP2IdxLocal));
                // and then we are also able to grab the corresponding vertices
                const auto &v1 = top->graph().vertices().at(v1Ix);
                const auto &v2 = top->graph().vertices().at(v2Ix);
                auto find1 = std::find(std::begin(v1.neighbors()), std::end(v1.neighbors()), v2Ix);
                auto find2 = std::find(std::begin(v2.neighbors()), std::end(v2.neighbors()), v1Ix);
                // making sure that the edge exists
                REQUIRE(find1 != std::end(v1.neighbors()));
                REQUIRE(find2 != std::end(v2.neighbors()));
            }

        }
    });

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
        for (time = 1; time < n_time_steps; ++time) {
            integrator->perform();
            topReactions->perform();
            forces->perform();
            reactions->perform();
            kernel->evaluateObservables(time);
            for (auto topPtr : kernel->stateModel().getTopologies()) {
                // check that all topologies are just containing T particles and their edges are also fine
                for (const auto &p : topPtr->fetchParticles()) {
                    REQUIRE(p.type() == ctx.particleTypes().idOf("T"));
                }
                for (auto edge : topPtr->graph().edges()) {
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
        top->addEdge(it.persistent_index(), it2.persistent_index());
        ++it;
        ++it2;
        top->addEdge(it.persistent_index(), it2.persistent_index());
    }

    simulation.addParticle("A", -2., 0, 0);
    simulation.addParticle("A", -3., 0, 0);
    simulation.addParticle("A", -4., 0, 0);
    simulation.addParticle("A", 2., 0, 0);
    simulation.addParticle("A", 3., 0, 0);
    simulation.addParticle("A", 4., 0, 0);

    simulation.run(6, 1.);

    const auto &type_registry = simulation.context().particleTypes();

    REQUIRE(simulation.context().topologyRegistry().isSpatialReactionType("A"));
    REQUIRE(simulation.context().topologyRegistry().isSpatialReactionType("end"));
    REQUIRE_FALSE(simulation.context().topologyRegistry().isSpatialReactionType("middle"));
    REQUIRE(simulation.context().calculateMaxCutoff() == 1.5);

    REQUIRE(simulation.currentTopologies().size() == 1);
    auto chainTop = simulation.currentTopologies().at(0);
    REQUIRE(chainTop->nParticles() == 3 /*original topology particles*/ + 6 /*attached particles*/);

    auto top_particles = simulation.stateModel().getParticlesForTopology(*chainTop);

    bool foundEndVertex{false};
    // check that graph is indeed linear
    for (std::size_t idx = 0; idx < chainTop->graph().vertices().size() && !foundEndVertex; ++idx) {
        auto prev_neighborIt = std::next(chainTop->graph().vertices().begin(), idx);
        auto prev_neighbor = prev_neighborIt.persistent_index();
        const auto &v_end = *prev_neighborIt;

        if (chainTop->particleForVertex(v_end).type() == type_registry.idOf("end")) {
            foundEndVertex = true;

            REQUIRE(v_end.neighbors().size() == 1);
            REQUIRE(top_particles.at(idx).type() == type_registry.idOf("end"));

            using flouble = readdy::fp::FloatingPoint<readdy::scalar>;
            flouble x_end(top_particles.at(idx).pos().x);
            flouble y_end(top_particles.at(idx).pos().y);
            flouble z_end(top_particles.at(idx).pos().z);
            // the end particle of our topology sausage should be either at x=4 or x=-4
            REQUIRE((x_end.AlmostEquals(flouble{4.}) || x_end.AlmostEquals(flouble{-4.})));
            REQUIRE(y_end.AlmostEquals(flouble{0.})); // no diffusion going on
            REQUIRE(z_end.AlmostEquals(flouble{0.})); // no diffusion going on

            auto factor = x_end.AlmostEquals(flouble{4.}) ? 1. : -1.;

            // walk along topology sausage, check end particles are always at +-4, the other ones are of type middle
            auto it = chainTop->graph().begin();
            auto endIt1 = it;
            auto endIt2 = it;
            bool firstEnd = true;
            for (; it != chainTop->graph().end(); ++it) {
                if (chainTop->particleForVertex(it.persistent_index()).type() == type_registry.idOf("end")) {
                    if (firstEnd) {
                        endIt1 = it;
                        firstEnd = false;
                    } else {
                        endIt2 = it;
                    }
                }
            }
            REQUIRE(endIt1 != chainTop->graph().end());
            REQUIRE(endIt2 != chainTop->graph().end());
            REQUIRE(chainTop->graph().graphDistance(endIt1, endIt2) == 8);

            if (flouble(chainTop->particleForVertex(endIt1.persistent_index()).pos().x).AlmostEquals(flouble(-4))) {
                REQUIRE(flouble(chainTop->particleForVertex(endIt2.persistent_index()).pos().x).AlmostEquals(
                        flouble(4)));
            } else if (flouble(chainTop->particleForVertex(endIt1.persistent_index()).pos().x).AlmostEquals(
                    flouble(4))) {
                REQUIRE(flouble(chainTop->particleForVertex(endIt2.persistent_index()).pos().x).AlmostEquals(
                        flouble(-4)));
                std::swap(endIt1, endIt2);
            } else {
                FAIL("the ends must be at {-4, 4}");
            }

            // now endIt1 is at -4 and endIt2 is at 4
            REQUIRE(chainTop->graph().vertices().at(endIt1.persistent_index()).neighbors().size() == 1);
            auto nneighbor = chainTop->graph().vertices().at(endIt1.persistent_index()).neighbors()[0];
            REQUIRE(chainTop->particleForVertex(nneighbor).type() == type_registry.idOf("middle"));
            REQUIRE(flouble(chainTop->particleForVertex(nneighbor).pos().x).AlmostEquals(flouble(-3)));
            REQUIRE(chainTop->graph().vertices().at(nneighbor).neighbors().size() == 2);

            auto neighborIx =
                    chainTop->graph().vertices().at(nneighbor).neighbors()[0] == endIt1.persistent_index() ? 1 : 0;
            auto prevNeighbor = nneighbor;
            nneighbor = chainTop->graph().vertices().at(nneighbor).neighbors()[neighborIx];
            for (int i = -2; i <= 3; ++i) {
                REQUIRE(chainTop->particleForVertex(nneighbor).type() == type_registry.idOf("middle"));
                REQUIRE(flouble(chainTop->particleForVertex(nneighbor).pos().x).AlmostEquals(flouble(i)));
                REQUIRE(chainTop->graph().vertices().at(nneighbor).neighbors().size() == 2);
                neighborIx = chainTop->graph().vertices().at(nneighbor).neighbors()[0] == prevNeighbor ? 1 : 0;
                prevNeighbor = nneighbor;
                nneighbor = chainTop->graph().vertices().at(nneighbor).neighbors()[neighborIx];
            }
            REQUIRE(chainTop->particleForVertex(nneighbor).type() == type_registry.idOf("end"));
            REQUIRE(flouble(chainTop->particleForVertex(nneighbor).pos().x).AlmostEquals(flouble(4)));
            REQUIRE(chainTop->graph().vertices().at(nneighbor).neighbors().size() == 1);
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
        ctx.potentials().addCylinder("A", 100., {0., 0., 0.}, {1., 0., 0.}, 0.01, true);
        ctx.potentials().addCylinder("head", 100., {0., 0., 0.}, {1., 0., 0.}, 0.01, true);
        ctx.potentials().addCylinder("tail", 100., {0., 0., 0.}, {1., 0., 0.}, 0.01, true);

        ctx.potentials().addBox("head", 10., {-4., -12.5, -12.5}, {0.000001, 25., 25.});

        std::vector<readdy::model::Particle> particles{
                {-4., 0., 0., types.idOf("head")},
                {-2., 0., 0., types.idOf("A")},
                {2.,  0., 0., types.idOf("A")},
                {4.,  0., 0., types.idOf("tail")},
        };
        REQUIRE(particles.size() == 4);

        auto graphTop = stateModel.addTopology(topReg.idOf("T1"), particles);
        {
            auto &graph = graphTop->graph();
            for (std::size_t i = 0; i < 3; ++i) {
                graphTop->addEdge({i}, {i + 1});
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
            for (std::size_t t = 1; t < nSteps + 1; t++) {
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
            for (std::size_t t = 1; t < nSteps + 1; t++) {
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

TEMPLATE_TEST_CASE("Helix grows by spatial topology reactions", "[!hide][integration]", SingleCPU, CPU) {
    GIVEN("An initial polymer B-A-A-...-A-C with helical structure") {
        std::vector<readdy::Vec3> initPos23 =
                {{-1.75373831, -2.14809599, 3.09407703},
                 {-0.75373831, -2.14809599, 3.09407703},
                 {0.13171772,  -1.68337282, 3.09407703},
                 {0.7030635,   -0.86664044, 3.0133787},
                 {0.8416891,   0.09866305,  2.79208227},
                 {0.54051056,  0.97080799,  2.40652464},
                 {-0.09420439, 1.53207468,  1.87537811},
                 {-0.87022492, 1.64486677,  1.25483795},
                 {-1.55890118, 1.2871743,   0.62413906},
                 {-1.95409477, 0.55824634,  0.06513454},
                 {-1.92530725, -0.34698867, -0.35880043},
                 {-1.45131376, -1.18816226, -0.61910247},
                 {-0.62563341, -1.74141669, -0.72938167},
                 {0.36757365,  -1.85710024, -0.74191387},
                 {1.30097023,  -1.49833734, -0.73416743},
                 {1.96263423,  -0.75053257, -0.7888372},
                 {2.21067781,  0.20064083,  -0.97253135},
                 {2.00981645,  1.11709499,  -1.3186043},
                 {1.44046291,  1.76969029,  -1.81856021},
                 {0.67800228,  1.99729206,  -2.42424266},
                 {-0.05241042, 1.74830016,  -3.06024607},
                 {-0.53387995, 1.09394804,  -3.64335016},
                 {-0.61367178, 0.20994353,  -4.10396854}};

        readdy::scalar angle = 2. * M_PI / 13;
        readdy::scalar dihedral = -10. * M_PI / 180.;

        readdy::model::Context ctx;
        ctx.boxSize() = {50., 50., 50.};
        ctx.particleTypes().add("A", 0.1, readdy::model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("B", 0.1, readdy::model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("C", 0.1, readdy::model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("S", 1000., readdy::model::particleflavor::TOPOLOGY);
        ctx.topologyRegistry().addType("helix");
        ctx.topologyRegistry().configureBondPotential("A", "A", {1000., 1.});
        ctx.topologyRegistry().configureBondPotential("B", "A", {1000., 1.});
        ctx.topologyRegistry().configureBondPotential("C", "A", {1000., 1.});
        ctx.topologyRegistry().configureAnglePotential("A", "A", "A", {1000., M_PI - angle});
        ctx.topologyRegistry().configureAnglePotential("B", "A", "A", {1000., M_PI - angle});
        ctx.topologyRegistry().configureAnglePotential("A", "A", "C", {1000., M_PI - angle});
        ctx.topologyRegistry().configureTorsionPotential("A", "A", "A", "A", {1000., 1., -dihedral});
        ctx.topologyRegistry().configureTorsionPotential("B", "A", "A", "A", {1000., 1., -dihedral});
        ctx.topologyRegistry().configureTorsionPotential("A", "A", "A", "C", {1000., 1., -dihedral});
        ctx.topologyRegistry().addSpatialReaction("Fusion: helix(C)+(S) -> helix(A--C)", 10000., 10.);

        readdy::Simulation simulation(create<TestType>(), ctx);
        WHEN("the helix elongates by attaching 20 segments from substrate particles") {
            {
                std::vector<readdy::model::Particle> topParticles;
                topParticles.reserve(23);
                topParticles.emplace_back(initPos23.front(), ctx.particleTypes().idOf("B"));
                for (std::size_t i = 1; i < initPos23.size() - 1; ++i) {
                    topParticles.emplace_back(initPos23[i], ctx.particleTypes().idOf("A"));
                }
                topParticles.emplace_back(initPos23.back(), ctx.particleTypes().idOf("C"));
                auto top = simulation.addTopology("helix", topParticles);
                auto &graph = top->graph();
                for (std::size_t i = 0; i < topParticles.size() - 1; ++i) {
                    top->addEdgeBetweenParticles(i, i + 1);
                }

                for (std::size_t i = 0; i < 20; ++i) {
                    readdy::scalar x = readdy::model::rnd::uniform_real(-25., 25.);
                    readdy::scalar y = readdy::model::rnd::uniform_real(-25., 25.);
                    readdy::scalar z = readdy::model::rnd::uniform_real(-25., 25.);
                    simulation.addParticle("S", x, y, z);
                }
            }

            std::vector<readdy::scalar> distances;
            auto callback = [&distances, &ctx](const std::vector<readdy::Vec3> &dist) {
                const auto d = readdy::bcs::dist(dist.at(0), dist.at(1), ctx.boxSize(),
                                                 ctx.periodicBoundaryConditions());
                distances.push_back(d);
            };
            auto obs = simulation.observe().positions(1000, {"B", "C"}, callback);
            simulation.registerObservable(std::move(obs));

            std::size_t nSteps = 2400000;
            readdy::scalar dt = 4e-5;
            simulation.run(nSteps, dt);

            REQUIRE(distances.size() > 0);

            THEN("the average end-to-end distance after absorption of all substrate particles is roughly ~14.4") {
                // extract recorded distances after a certain number of steps
                // when all substrates are most likely absorbed
                auto n = distances.size();
                auto first = distances.begin() + 2 * n / 10; // i.e. after 20% of the time
                std::vector<readdy::scalar> distances2(first, distances.end());
                auto n2 = distances2.size();
                auto mean = std::accumulate(distances2.begin(), distances2.end(), 0.);
                mean /= static_cast<readdy::scalar>(n2);
                readdy::scalar expected = 14.434;
                auto relativeErr = std::abs(mean - expected) / expected;
                // allow for 5% deviation. Usually the observed mean is expected to deviate about ~0.5%
                CHECK(relativeErr < 0.05);
            }
        }
    }
}

TEMPLATE_TEST_CASE("Helix grows by structural topology reactions", "[!hide][integration]", SingleCPU, CPU) {
    GIVEN("An initial polymer B-A-A-...-A-C with helical structure") {
        std::vector<readdy::Vec3> initPos23 = {{-1.75373831, -2.14809599, 3.09407703},
                                       {-0.75373831, -2.14809599, 3.09407703},
                                       {0.13171772,  -1.68337282, 3.09407703},
                                       {0.7030635,   -0.86664044, 3.0133787},
                                       {0.8416891,   0.09866305,  2.79208227},
                                       {0.54051056,  0.97080799,  2.40652464},
                                       {-0.09420439, 1.53207468,  1.87537811},
                                       {-0.87022492, 1.64486677,  1.25483795},
                                       {-1.55890118, 1.2871743,   0.62413906},
                                       {-1.95409477, 0.55824634,  0.06513454},
                                       {-1.92530725, -0.34698867, -0.35880043},
                                       {-1.45131376, -1.18816226, -0.61910247},
                                       {-0.62563341, -1.74141669, -0.72938167},
                                       {0.36757365,  -1.85710024, -0.74191387},
                                       {1.30097023,  -1.49833734, -0.73416743},
                                       {1.96263423,  -0.75053257, -0.7888372},
                                       {2.21067781,  0.20064083,  -0.97253135},
                                       {2.00981645,  1.11709499,  -1.3186043},
                                       {1.44046291,  1.76969029,  -1.81856021},
                                       {0.67800228,  1.99729206,  -2.42424266},
                                       {-0.05241042, 1.74830016,  -3.06024607},
                                       {-0.53387995, 1.09394804,  -3.64335016},
                                       {-0.61367178, 0.20994353,  -4.10396854}};
        readdy::scalar angle = 2. * M_PI / 13;
        readdy::scalar dihedral = -10. * M_PI / 180.;
        readdy::model::Context ctx;
        ctx.boxSize() = {50., 50., 50.};
        ctx.particleTypes().add("A", 0.1, readdy::model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("B", 0.1, readdy::model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("C", 0.1, readdy::model::particleflavor::TOPOLOGY);
        ctx.topologyRegistry().addType("helix");
        ctx.topologyRegistry().configureBondPotential("A", "A", {1000., 1.});
        ctx.topologyRegistry().configureBondPotential("B", "A", {1000., 1.});
        ctx.topologyRegistry().configureBondPotential("C", "A", {1000., 1.});
        ctx.topologyRegistry().configureAnglePotential("A", "A", "A", {1000., M_PI - angle});
        ctx.topologyRegistry().configureAnglePotential("B", "A", "A", {1000., M_PI - angle});
        ctx.topologyRegistry().configureAnglePotential("A", "A", "C", {1000., M_PI - angle});
        ctx.topologyRegistry().configureTorsionPotential("A", "A", "A", "A", {1000., 1., -dihedral});
        ctx.topologyRegistry().configureTorsionPotential("B", "A", "A", "A", {1000., 1., -dihedral});
        ctx.topologyRegistry().configureTorsionPotential("A", "A", "A", "C", {1000., 1.,-dihedral});
        namespace rmt = readdy::model::top;
        auto reactionFunction = [&ctx](rmt::GraphTopology &topology) -> rmt::reactions::Recipe {
            rmt::reactions::Recipe recipe(topology);
            auto &graph = topology.graph();
            for (auto vIt = graph.begin(); vIt != graph.end(); ++vIt) {
                auto p = topology.particleForVertex(vIt.persistent_index());
                if (p.type() == ctx.particleTypes().idOf("C")) {
                    readdy::Vec3 pos1;
                    for (auto &neighbor : vIt->neighbors()) { pos1 = topology.particleForVertex(neighbor).pos(); }
                    auto pos2 = topology.particleForVertex(vIt.persistent_index()).pos();
                    auto newPosition = pos2 + (pos2 - pos1);
                    auto ida = ctx.particleTypes().idOf("A");
                    recipe.changeParticleType(vIt.persistent_index(), ida);
                    std::string ctype = "C";
                    recipe.appendNewParticle({vIt.persistent_index()}, ctype, newPosition);
                }
            }
            return recipe;
        };
        auto rateFunction = [](const rmt::GraphTopology &topology) -> readdy::scalar {
            auto nVertices = topology.graph().vertices().size();
            if (nVertices < 43) { return 1e2; } else { return 0.; }
        };
        readdy::model::top::reactions::StructuralTopologyReaction reaction("append", reactionFunction, rateFunction);
        std::string type = "helix";
        ctx.topologyRegistry().addStructuralReaction(type, reaction);
        readdy::Simulation simulation(create<TestType>(), ctx);
        WHEN("the helix elongates by attaching 20 segments from substrate particles") {
            {
                std::vector<readdy::model::Particle> topParticles;
                topParticles.reserve(23);
                topParticles.emplace_back(initPos23.front(), ctx.particleTypes().idOf("B"));
                for (std::size_t i = 1; i < initPos23.size() - 1; ++i) {
                    topParticles.emplace_back(initPos23[i], ctx.particleTypes().idOf("A"));
                }
                topParticles.emplace_back(initPos23.back(), ctx.particleTypes().idOf("C"));
                auto top = simulation.addTopology("helix", topParticles);
                auto &graph = top->graph();
                for (std::size_t i = 0; i < topParticles.size() - 1; ++i) {
                    top->addEdgeBetweenParticles(i, i + 1);
                }
            }
            std::vector<readdy::scalar> distances;
            auto callback = [&distances, &ctx](const std::vector<readdy::Vec3> &dist) {
                const auto d = readdy::bcs::dist(dist.at(0), dist.at(1), ctx.boxSize(), ctx.periodicBoundaryConditions());
                distances.push_back(d);
            };
            auto obs = simulation.observe().positions(1000, {"B", "C"}, callback);

            simulation.registerObservable(std::move(obs));
            std::size_t nSteps = 2400000;
            readdy::scalar dt = 4e-5;
            simulation.run(nSteps, dt);
            REQUIRE(!distances.empty());
            THEN("the average end-to-end distance after absorption of all substrate particles is roughly ~14.4") {
                /* extract recorded distances after a certain number of steps when all substrates are most likely absorbed */
                auto n = distances.size();
                auto first = distances.begin() + 2 * n / 10; /* i.e. after 20% of the time*/
                std::vector<readdy::scalar> distances2(first, distances.end());
                auto n2 = distances2.size();
                auto mean = std::accumulate(distances2.begin(), distances2.end(), 0.);
                mean /= static_cast<readdy::scalar>(n2);
                readdy::scalar expected = 14.434;
                auto relativeErr = std::abs(mean - expected) / expected;
                /* allow for 5% deviation. Usually the observed mean is expected to deviate about ~0.5% */
                CHECK(relativeErr < 0.05);
            }
        }
    }
}

TEMPLATE_TEST_CASE("Particles form complexes with predetermined number of bonds", "[!hide][integration]", SingleCPU, CPU) {
    /* Credits go to Moritz FP Becker for this test case. */
    readdy::model::Context ctx;
    ctx.boxSize() = {25, 25, 25};
    ctx.topologyRegistry().addType("Complex");
    ctx.particleTypes().add("A", 1., readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("Aa", 1., readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("Aaa", 1., readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("Aaaa", 1., readdy::model::particleflavor::TOPOLOGY);

    ctx.topologyRegistry().configureBondPotential("Aa", "Aa", {.forceConstant = 0., .length = 1.});
    ctx.topologyRegistry().configureBondPotential("Aa", "Aaa", {.forceConstant = 0., .length = 1.});
    ctx.topologyRegistry().configureBondPotential("Aa", "Aaaa", {.forceConstant = 0., .length = 1.});
    ctx.topologyRegistry().configureBondPotential("Aaa", "Aaa", {.forceConstant = 0., .length = 1.});
    ctx.topologyRegistry().configureBondPotential("Aaa", "Aaaa", {.forceConstant = 0., .length = 1.});
    ctx.topologyRegistry().configureBondPotential("Aaaa", "Aaaa", {.forceConstant = 0., .length = 1.});

    ctx.topologyRegistry().addSpatialReaction("bind_1: Complex(A) + Complex(A) -> Complex(Aa--Aa) [self=true]", 1., 1.);
    ctx.topologyRegistry().addSpatialReaction("bind_2: Complex(A) + Complex(Aa) -> Complex(Aa--Aaa) [self=true]", 1., 1.);
    ctx.topologyRegistry().addSpatialReaction("bind_3: Complex(A) + Complex(Aaa) -> Complex(Aa--Aaaa) [self=true]", 1., 1.);
    ctx.topologyRegistry().addSpatialReaction("bind_4: Complex(Aa) + Complex(Aa) -> Complex(Aaa--Aaa) [self=true]", 1., 1.);
    ctx.topologyRegistry().addSpatialReaction("bind_5: Complex(Aa) + Complex(Aaa) -> Complex(Aaa--Aaaa) [self=true]", 1., 1.);
    ctx.topologyRegistry().addSpatialReaction("bind_6: Complex(Aaa) + Complex(Aaa) -> Complex(Aaaa--Aaaa) [self=true]", 1., 1.);

    readdy::Simulation simulation(create<TestType>(), ctx);

    auto nParticles = 1000;
    for (auto i = 0U; i < nParticles; ++i) {
        readdy::Vec3 pos {readdy::model::rnd::uniform_real(-25./2., 25./2.),
                          readdy::model::rnd::uniform_real(-25./2., 25./2.),
                          readdy::model::rnd::uniform_real(-25./2., 25./2.)};
        readdy::model::Particle A {pos, ctx.particleTypes().idOf("A")};
        simulation.addTopology("Complex", {A});
    }
    simulation.run(10000, 1e-2);
    auto topologies = simulation.currentTopologies();

    std::map<std::size_t, std::tuple<std::string, std::size_t, std::size_t>> edgeCounts;
    std::size_t tix = 0;
    for(auto top : topologies) {
        for(const auto & edge : top->graph().edges()) {
            const auto &[pix1, pix2] = edge;

            {
                auto p1 = top->particleForVertex(pix1);
                auto c1 = std::get<1>(edgeCounts[p1.id()]);
                edgeCounts[p1.id()] = std::make_tuple(ctx.particleTypes().nameOf(p1.type()), c1 + 1, tix);
            }
            {
                auto p2 = top->particleForVertex(pix2);
                auto c2 = std::get<1>(edgeCounts[p2.id()]);
                edgeCounts[p2.id()] = std::make_tuple(ctx.particleTypes().nameOf(p2.type()), c2 + 1, tix);
            }
        }
        ++tix;
    }

    auto particles = simulation.stateModel().getParticles();
    REQUIRE(particles.size() == nParticles);
    int singleParticles = 0;
    for (const auto& p : particles) {
        if(edgeCounts.find(p.id()) == edgeCounts.end()) {
            REQUIRE(ctx.particleTypes().nameOf(p.type()) == "A");
            ++singleParticles;
        }
    }

    for(auto top : simulation.currentTopologies()) {
        for(auto v : top->graph().vertices()) {
            auto p = top->particleForVertex(v);
            REQUIRE(v.neighbors().size() == ctx.particleTypes().nameOf(p.type()).size() - 1);
        }
    }

    REQUIRE(edgeCounts.size() + singleParticles == nParticles);
    for (auto entry : edgeCounts) {
        REQUIRE(std::get<0>(entry.second).size() - 1 == std::get<1>(entry.second));
    }
}
