/**
 * << detailed description >>
 *
 * @file SingleCPUTestParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.06.16
 */

#include <gtest/gtest.h>
#include <readdy/common/make_unique.h>
#include <readdy/kernel/singlecpu/model/SingleCPUParticleData.h>
#include <boost/log/trivial.hpp>


using namespace readdy::kernel::singlecpu::model;

TEST(ParticleData, initialization) {
    // everything should be deactivated
    auto data = std::make_unique<SingleCPUParticleData>(10);
    auto&& it = data->getDeactivatedParticles()->begin();

    for (size_t i = 0; i < 10; i++) {
        EXPECT_EQ(*it, i);
        ++it;
    }
}

TEST(ParticleData, additionOfParticles) {
    // initial capacity of 3
    auto data = std::make_unique<SingleCPUParticleData>(3);
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 3);

    data->addParticle({1, 1, 1, 5});
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 2);

    data->addParticle({2, 2, 2, 6});
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 1);

    data->addParticle({3, 3, 3, 7});
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 0);

    data->addParticle({4, 4, 4, 8});
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 0);

    int i = 0;
    for (auto &&it = data->begin_positions(); it != data->end_positions(); ++it) {
        i++;
        EXPECT_TRUE(*it == readdy::model::Vec3(i, i, i));
    }

    i = 5;
    for (auto &&it = data->begin_types(); it != data->end_types(); ++it) {
        EXPECT_TRUE((*it) == i);
        i++;
    }
}

TEST(ParticleData, removalOfParticles) {
    auto data = std::make_unique<SingleCPUParticleData>(5);
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 5);
    for (auto &&i = 0; i < 100; i++) {
        data->addParticle({(double) i, (double) i, (double) i, 5});
    }
    EXPECT_TRUE(data->size() == 100);
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 0);
    int i = 0;
    for (auto &&it = data->begin_positions(); it != data->end_positions(); ++it) {
        EXPECT_TRUE(*it == readdy::model::Vec3(i, i, i));
        i++;
    }

    data->removeParticle(2);
    EXPECT_TRUE(data->size() == 99);
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 1);
    EXPECT_TRUE(*data->getDeactivatedParticles()->begin() == 2);
    // we expect an offset by 1 for every particle that comes after the 3rd
    EXPECT_TRUE(*(data->begin_positions() + 10) == readdy::model::Vec3(11, 11, 11));
    EXPECT_TRUE(*(data->begin_positions() + 1) == readdy::model::Vec3(1, 1, 1));

    // take up the space again
    data->addParticle({.5, .5, .5, 1});
    EXPECT_TRUE(data->size() == 100);
    EXPECT_TRUE(data->getDeactivatedParticles()->size() == 0);
    EXPECT_TRUE(*(data->begin_positions() + 2) == readdy::model::Vec3(.5, .5, .5));
}

TEST(ParticleData, empty) {
    auto data = std::make_unique<SingleCPUParticleData>(3);
    EXPECT_TRUE(data->empty());
    data->addParticle({4, 4, 4, 4});
    EXPECT_FALSE(data->empty());
}