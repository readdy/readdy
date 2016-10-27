/**
 * << detailed description >>
 *
 * @file TestVec3.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <gtest/gtest.h>
#include <readdy/model/Vec3.h>
#include <readdy/common/logging.h>

namespace {

using vec_t = readdy::model::Vec3;
namespace log = readdy::log;

    TEST(Vec3, SizeOfVec3) {
        vec_t vec(0,0,0);
        EXPECT_EQ(8*3, sizeof(vec));
    }
}