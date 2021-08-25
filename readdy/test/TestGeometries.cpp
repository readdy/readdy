//
// Created by mho on 8/25/21.
//

#include <catch2/catch.hpp>
#include <readdy/model/geometry.h>

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
}

TEST_CASE("Sphere geometry", "[geometries]") {
    model::geometry::Sphere<scalar> sphere {
        .center = {1., 1., 1.},
        .radius = 5.
    };
    REQUIRE(sphere.contains<true>({0., 0., 0.}));
    REQUIRE(!sphere.contains<false>({0., 0., 0.}));

    REQUIRE(sphere.contains<true>({1., 1., 1.}));
    REQUIRE(!sphere.contains<false>({1., 1., 1.}));

    REQUIRE(!sphere.contains<true>({20., 20., 20.}));
    REQUIRE(sphere.contains<false>({20., 20., 20.}));
}
