/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * @file TestMatrix33.cpp
 * @author clonker
 * @date 1/17/18
 */

#include <catch2/catch.hpp>

#include <readdy/common/numeric.h>
#include <readdy/common/common.h>

using mat = readdy::Matrix33;
using vec = readdy::Vec3;

SCENARIO("Working with matrices") {
    GIVEN("Two vectors") {
        vec u(1, 2, 3);
        vec v(3, 2, 1);
        WHEN("taking the outer product") {
            auto outer = readdy::math::outerProduct<mat>(u, v);
            THEN("a matrix is yielded") {
                REQUIRE(outer.at(0, 0) == u[0] * v[0]);
                REQUIRE(outer.at(0, 1) == u[0] * v[1]);
                REQUIRE(outer.at(0, 2) == u[0] * v[2]);
                REQUIRE(outer.at(1, 0) == u[1] * v[0]);
                REQUIRE(outer.at(1, 1) == u[1] * v[1]);
                REQUIRE(outer.at(1, 2) == u[1] * v[2]);
                REQUIRE(outer.at(2, 0) == u[2] * v[0]);
                REQUIRE(outer.at(2, 1) == u[2] * v[1]);
                REQUIRE(outer.at(2, 2) == u[2] * v[2]);
            }
        }
    }
    GIVEN("A 3x3 matrix with some elements") {
        mat m{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
        THEN("the elements should be retrievable") {
            REQUIRE(m.at(0, 0) == 1);
            REQUIRE(m.at(0, 1) == 2);
            REQUIRE(m.at(0, 2) == 3);
            REQUIRE(m.at(1, 0) == 4);
            REQUIRE(m.at(1, 1) == 5);
            REQUIRE(m.at(1, 2) == 6);
            REQUIRE(m.at(2, 0) == 7);
            REQUIRE(m.at(2, 1) == 8);
            REQUIRE(m.at(2, 2) == 9);
            REQUIRE(m.at(0, 0) == m.data()[0]);
            REQUIRE(m.at(0, 1) == m.data()[1]);
            REQUIRE(m.at(0, 2) == m.data()[2]);
            REQUIRE(m.at(1, 0) == m.data()[3]);
            REQUIRE(m.at(1, 1) == m.data()[4]);
            REQUIRE(m.at(1, 2) == m.data()[5]);
            REQUIRE(m.at(2, 0) == m.data()[6]);
            REQUIRE(m.at(2, 1) == m.data()[7]);
            REQUIRE(m.at(2, 2) == m.data()[8]);
        }
        THEN("the dims are 3x3") {
            REQUIRE(mat::n() == 3);
            REQUIRE(mat::m() == 3);
        }
        WHEN("taking a copy") {
            auto copy = m;
            mat copy2(m);
            THEN("modifying the copy should have no influence on the original matrix") {
                copy.at(0, 1) = 5;
                copy2.at(0, 1) = 6;

                REQUIRE(m.at(0, 1) == 2);
                REQUIRE(copy.at(0, 1) == 5);
                REQUIRE(copy2.at(0, 1) == 6);
            }
        }
        WHEN("adding another matrix") {
            mat m2{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
            m += m2;
            auto m3 = m2 + m2;

            THEN("the result should be the sum of the two original matrices") {
                CHECK(m.at(0, 0) == 2);
                CHECK(m.at(0, 1) == 4);
                CHECK(m.at(0, 2) == 6);
                CHECK(m.at(1, 0) == 8);
                CHECK(m.at(1, 1) == 10);
                CHECK(m.at(1, 2) == 12);
                CHECK(m.at(2, 0) == 14);
                CHECK(m.at(2, 1) == 16);
                CHECK(m.at(2, 2) == 18);
                CHECK(m3.at(0, 0) == 2);
                CHECK(m3.at(0, 1) == 4);
                CHECK(m3.at(0, 2) == 6);
                CHECK(m3.at(1, 0) == 8);
                CHECK(m3.at(1, 1) == 10);
                CHECK(m3.at(1, 2) == 12);
                CHECK(m3.at(2, 0) == 14);
                CHECK(m3.at(2, 1) == 16);
                CHECK(m3.at(2, 2) == 18);
                CHECK(m == m3);
            }
        }
        WHEN("scaling the matrix") {
            auto m2 = m * 2;
            auto m3 = 2 * m;
            m *= 2;

            THEN("the result should be the original matrix but scaled") {
                REQUIRE(m.at(0, 0) == 2);
                REQUIRE(m.at(0, 1) == 4);
                REQUIRE(m.at(0, 2) == 6);
                REQUIRE(m.at(1, 0) == 8);
                REQUIRE(m.at(1, 1) == 10);
                REQUIRE(m.at(1, 2) == 12);
                REQUIRE(m.at(2, 0) == 14);
                REQUIRE(m.at(2, 1) == 16);
                REQUIRE(m.at(2, 2) == 18);
                REQUIRE(m == m2);
                REQUIRE(m == m3);
            }
        }

        WHEN("given two other matrices") {
            mat m2{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
            mat m3{{{2, 2, 3, 4, 5, 6, 7, 8, 9}}};

            THEN("equality can be checked") {
                REQUIRE(m == m);
                REQUIRE(m == m2);
                REQUIRE_FALSE(m == m3);
            }
        }
    }
}
