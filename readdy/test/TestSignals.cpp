//
// Created by mho on 17/10/2016.
//

#include <gtest/gtest.h>
#include <readdy/common/signals.h>

namespace sig = readdy::signals;

namespace {

struct TestSlot {

    TestSlot(const std::function<void()> &callback) : callback(callback) {}

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