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


//
// Created by mho on 17/10/2016.
//

#include <gtest/gtest.h>
#include <readdy/common/signals.h>

namespace sig = readdy::signals;

namespace {

struct TestSlot {

    explicit TestSlot(const std::function<void()> &callback) : callback(callback) {}

    void operator()() {
        callback();
    }

private:
    std::function<void()> callback;
};

TEST(TestSignals, TestSimpleSignalConnectDisconnect) {
    sig::signal<void(int)> signal;

    {
        auto connection = signal.connect([](int i) {});
        EXPECT_EQ(signal.n_slots(), 1);
        connection.disconnect();
        EXPECT_EQ(signal.n_slots(), 0);
    }
    {
        auto scoped = signal.connect_scoped([](int i) {});
        EXPECT_EQ(signal.n_slots(), 1);
    }
    EXPECT_EQ(signal.n_slots(), 0);
}

TEST(TestSignals, TestStructSignalConnectDisconnect) {
    sig::signal<void()> signal;
    unsigned int called = 0;
    {
        TestSlot slot{[&called] {
            ++called;
        }};
        {
            auto scoped = signal.connect_scoped(slot);
            EXPECT_EQ(signal.n_slots(), 1);
            signal.fire_signal();
        }
        slot();
        EXPECT_EQ(signal.n_slots(), 0);
        EXPECT_EQ(called, 2);
    }
}

}