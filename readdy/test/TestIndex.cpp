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
 * << detailed description >>
 *
 * @file TextIndex.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.09.17
 * @copyright BSD-3
 */

#include <catch2/catch.hpp>

#include <readdy/common/common.h>
#include <readdy/common/Index.h>

using namespace readdy;

TEST_CASE("Test the multidim-index class.", "[index]") {
    SECTION("1D") {
        util::Index1D index (std::array<std::size_t, 1>{5_z});
        CHECK(index[0] == 5);
        CHECK(index.get<0>() == 5);
        CHECK(index.nElements() == 5);
        CHECK(index(3_z) == 3);
        CHECK(index.inverse(3)[0] == 3);
    }
    SECTION("2D") {
        auto nrows = 5_z;
        auto ncols = 6_z;
        util::Index2D index (std::array<std::size_t, 2>{nrows, ncols});
        CHECK(index[0] == nrows);
        CHECK(index[1] == ncols);
        CHECK(index.get<0>() == nrows);
        CHECK(index.get<1>() == ncols);
        CHECK(index.nElements() == nrows * ncols);
        CHECK(index(0, 3) == 3);
        CHECK(index(1, 3) == 3 + ncols);
        CHECK(index(4, 5) == 5 + 4 * ncols);
        CHECK(index.inverse(3)[0] == 0);
        CHECK(index.inverse(3)[1] == 3);
        CHECK(index.inverse(3+ncols)[0] == 1);
        CHECK(index.inverse(3+ncols)[1] == 3);
        CHECK(index.inverse(5+4*ncols)[0] == 4);
        CHECK(index.inverse(5+4*ncols)[1] == 5);
    }
    SECTION("3D") {
        auto width = 6_z;
        auto height = 5_z;
        auto depth = 6_z;
        util::Index3D index (std::array<std::size_t, 3>{width, height, depth});
        CHECK(index[0] == index.get<0>());
        CHECK(index[1] == index.get<1>());
        CHECK(index[2] == index.get<2>());
        CHECK(index[0] == width);
        CHECK(index[1] == height);
        CHECK(index[2] == depth);
        CHECK(index.nElements() == width * height * depth);
        CHECK(index(0, 0, 3) == 3);
        CHECK(index(1, 2, 3) == 3 + depth * (2 + height * 1));
        util::Index3D::GridDims expected {{1_z,2_z,3_z}};
        CHECK(index.inverse(3 + depth * (2 + height * 1)) == expected);
    }
    SECTION("3D extensive") {
        auto width = 7_z;
        auto height = 13_z;
        auto depth = 5_z;

        util::Index3D index (std::array<std::size_t, 3>{width, height, depth});

        int n = 0;
        for(int i = 0; i < width; ++i) {
            for(int j = 0; j < height; ++j) {
                for(int k = 0; k < depth; ++k) {
                    auto x = index(i, j, k);
                    CHECK(x == n);
                    util::Index3D::GridDims expectedInverse {
                            {static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k)}
                    };
                    CHECK(index.inverse(x) == expectedInverse);
                    ++n;
                }
            }
        }
    }
    SECTION("5D") {
        using namespace readdy;
        auto d1 = 7_z;
        auto d2 = 13_z;
        auto d3 = 5_z;
        auto d4 = 11_z;
        auto d5 = 23_z;

        util::Index<5> index (std::array<std::size_t, 5>{d1, d2, d3, d4, d5});

        int n = 0;
        for(int i1 = 0; i1 < d1; ++i1) {
            for(int i2 = 0; i2 < d2; ++i2) {
                for(int i3 = 0; i3 < d3; ++i3) {
                    for(int i4 = 0; i4 < d4; ++i4) {
                        for(int i5 = 0; i5 < d5; ++i5) {
                            auto x = index(i1, i2, i3, i4, i5);
                            CHECK(x == n);
                            util::Index<5>::GridDims expectedInverse {
                                    {
                                            static_cast<std::size_t>(i1), static_cast<std::size_t>(i2),
                                            static_cast<std::size_t>(i3), static_cast<std::size_t>(i4),
                                            static_cast<std::size_t>(i5)
                                    }
                            };
                            CHECK(index.inverse(x) == expectedInverse);
                            ++n;
                        }
                    }
                }
            }
        }
    }
}
