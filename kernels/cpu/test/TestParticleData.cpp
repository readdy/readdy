/**
 * << detailed description >>
 *
 * @file TestFoo.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <readdy/common/logging.h>
#include <readdy/kernel/cpu/model/ParticleData.h>
#include "gtest/gtest.h"

namespace {
    TEST(TestParticleData, TestEntryBytesize) {
        readdy::kernel::cpu::model::ParticleData::Entry entry {};
        EXPECT_EQ(64, sizeof(entry)) << "an entry should have exactly 64 bytes";
    }
}