//
// Created by mho on 8/25/21.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <readdy/model/geometry.h>
#include <readdy/model/Context.h>
#include <readdy/api/Simulation.h>
#include <readdy/plugin/KernelProvider.h>

using namespace readdy;

TEST_CASE("Box geometry", "[geometry]") {
    model::geometry::Box<scalar> box {
        .v0 = {-1, -1, -1},
        .v1 = {1, 1, 1}
    };
    REQUIRE(box.contains<true>({0, 0, 0}));
    REQUIRE(!box.contains<false>({0, 0, 0}));

    for(auto notContainedVec : {Vec3{55, 0, 0}, Vec3{0, 55, 0}, Vec3{0, 0, 55}, Vec3{55, 55, 0}, Vec3{55, 0, 55}, Vec3{0, 55, 55}, Vec3{55, 55, 55}}) {
        REQUIRE(box.contains<false>(notContainedVec));
        REQUIRE(!box.contains<true>(notContainedVec));
        REQUIRE(box.smallestDifference<false>(notContainedVec) == Vec3{0, 0, 0});
    }
    REQUIRE(box.smallestDifference<true>({55, 0, 0}) == Vec3{54, 0, 0});
    REQUIRE(box.smallestDifference<true>({0, 55, 0}) == Vec3{0, 54, 0});
    REQUIRE(box.smallestDifference<true>({0, 0, 55}) == Vec3{0, 0, 54});

    REQUIRE(box.smallestDifference<false>({.1, 0, 0}) == Vec3{-.9, 0, 0});
    REQUIRE(box.smallestDifference<false>({0, .1, 0}) == Vec3{0, -.9, 0});
    REQUIRE(box.smallestDifference<false>({0, 0, .1}) == Vec3{0, 0, -.9});
    REQUIRE(box.smallestDifference<false>({.2, .1, 0}) == Vec3{-.8, 0, 0});
    REQUIRE(box.smallestDifference<false>({.2, .1, -.5}) == Vec3{0, 0, .5});

    readdy::model::Context ctx;
    ctx.boxSize() = {{10, 10, 10}};
    ctx.particleTypes().add("A", 1.);
    ctx.potentials().addHarmonicGeometry("A", 100., model::geometry::Box<scalar>{
        .v0 = {-3, -3, -3},
        .v1 = {3, 3, 3}
    }, false);
    readdy::Simulation sim {readdy::plugin::KernelProvider::getInstance().create("SingleCPU"), ctx};
    for(std::size_t i = 0; i < 5000; ++i) {
        sim.addParticle("A", 0, 0, 0);
    }
    sim.run(1000, 1e-3);
    auto positions = sim.getAllParticlePositions();
    for (auto pos : positions) {
        REQUIRE(box.contains<false>(pos));
    }
}

TEST_CASE("Sphere geometry", "[geometries]") {
    model::geometry::Sphere<scalar> sphere{
            .center = {1., 1., 1.},
            .radius = 5.
    };
    REQUIRE(sphere.contains<true>({0., 0., 0.}));
    REQUIRE(!sphere.contains<false>({0., 0., 0.}));

    REQUIRE(sphere.contains<true>({1., 1., 1.}));
    REQUIRE(!sphere.contains<false>({1., 1., 1.}));

    REQUIRE(!sphere.contains<true>({20., 20., 20.}));
    REQUIRE(sphere.contains<false>({20., 20., 20.}));


    REQUIRE(sphere.smallestDifference<true>({0, 0, 0}) == Vec3{0, 0, 0});
    REQUIRE(sphere.smallestDifference<false>({0, 0, 0}) == (std::sqrt(3) - 5) * Vec3{-1, -1, -1} / std::sqrt(3));
    auto delta = sphere.smallestDifference<false>({0, 0, 0});
    REQUIRE((Vec3{0, 0, 0} - delta - sphere.center).norm() == Catch::Approx(sphere.radius));
    REQUIRE((Vec3{2, -1, .3} - sphere.smallestDifference<false>({2, -1, .3}) - sphere.center).norm() == Catch::Approx(sphere.radius));

    REQUIRE(sphere.smallestDifference<true>({7, 1, 1}) == Vec3{1, 0, 0});
    REQUIRE(sphere.smallestDifference<true>({1, 7, 1}) == Vec3{0, 1, 0});
    REQUIRE(sphere.smallestDifference<true>({1, 1, 7}) == Vec3{0, 0, 1});

    readdy::model::Context ctx;
    ctx.boxSize() = {{13, 13, 13}};
    ctx.particleTypes().add("A", 1.);
    ctx.potentials().addHarmonicGeometry("A", 100., model::geometry::Sphere<scalar>{
            .center = sphere.center,
            .radius = 7
    }, false);
    readdy::Simulation sim {readdy::plugin::KernelProvider::getInstance().create("SingleCPU"), ctx};
    for(std::size_t i = 0; i < 5000; ++i) {
        sim.addParticle("A", 0, 0, 0);
    }
    sim.run(1000, 1e-3);
    auto positions = sim.getAllParticlePositions();
    for (auto pos : positions) {
        REQUIRE(sphere.contains<false>(pos));
    }
}
