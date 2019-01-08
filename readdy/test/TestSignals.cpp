#include <utility>

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


//
// Created by mho on 17/10/2016.
//


#include <catch2/catch.hpp>
#include <readdy/common/signals.h>

namespace sig = readdy::signals;

struct TestSlot {

    explicit TestSlot(std::function<void()> callback) : callback(std::move(callback)) {}

    void operator()() {
        callback();
    }

private:
    std::function<void()> callback;
};

TEST_CASE("Test the signals class") {
    SECTION("Connect and disconnect") {
        sig::signal<void(int)> signal;
        {
            auto connection = signal.connect([](int i) {});
            REQUIRE(signal.n_slots() == 1);
            connection.disconnect();
            REQUIRE(signal.n_slots() == 0);
        }
        {
            auto scoped = signal.connect_scoped([](int i) {});
            REQUIRE(signal.n_slots() == 1);
        }
        REQUIRE(signal.n_slots() == 0);
    }
    SECTION("Connect and disconnect with struct") {
        sig::signal<void()> signal;
        unsigned int called = 0;
        {
            TestSlot slot{[&called] {
                ++called;
            }};
            {
                auto scoped = signal.connect_scoped(slot);
                REQUIRE(signal.n_slots() == 1);
                signal.fire_signal();
            }
            slot();
            REQUIRE(signal.n_slots() == 0);
            REQUIRE(called == 2);
        }
    }
}
