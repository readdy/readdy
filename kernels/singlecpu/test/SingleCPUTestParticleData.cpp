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

    for (size_t i = 0; i < 10; i++) {
        EXPECT_TRUE(data->isMarkedForDeactivation(i));
    }
}

TEST(ParticleData, additionOfParticles) {
    // initial capacity of 3
    auto data = std::make_unique<SingleCPUParticleData>(3);
    EXPECT_TRUE(data->getNDeactivated() == 3);

    data->addParticle({1, 1, 1, 5});
    EXPECT_TRUE(data->getNDeactivated() == 2);

    data->addParticle({2, 2, 2, 6});
    EXPECT_TRUE(data->getNDeactivated() == 1);

    data->addParticle({3, 3, 3, 7});
    EXPECT_TRUE(data->getNDeactivated() == 0);

    data->addParticle({4, 4, 4, 8});
    EXPECT_TRUE(data->getNDeactivated() == 0);

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
    unsigned int n_particles = 15;
    auto data = std::make_unique<SingleCPUParticleData>(5);
    EXPECT_TRUE(data->getNDeactivated() == 5);
    for (auto &&i = 0; i < n_particles; i++) {
        data->addParticle({(double) i, (double) i, (double) i, 5});
    }
    EXPECT_TRUE(data->size() == n_particles);
    EXPECT_TRUE(data->getNDeactivated() == 0);
    int i = 0;
    for (auto &&it = data->begin_positions(); it != data->end_positions(); ++it) {
        EXPECT_TRUE(*it == readdy::model::Vec3(i, i, i));
        i++;
    }

    // 2nd particle -> n_particles-1th
    data->removeParticle(2);
    // 2nd particle -> n_particles-2nd
    data->removeParticle(2);
    // 2nd particle -> n_particles-3rd
    data->removeParticle(5);
    EXPECT_EQ(data->size(), n_particles - 3);
    EXPECT_EQ(data->getNDeactivated(), 3);

    EXPECT_EQ(*(data->begin_positions() + 10), readdy::model::Vec3(10, 10, 10));
    EXPECT_EQ(*(data->begin_positions() + 1), readdy::model::Vec3(1, 1, 1));
    EXPECT_EQ(*(data->begin_positions() + 2), readdy::model::Vec3(n_particles - 2, n_particles - 2, n_particles - 2));
    EXPECT_EQ(*(data->begin_positions() + 5), readdy::model::Vec3(n_particles - 3, n_particles - 3, n_particles - 3));

    // take up the space again
    data->addParticle({.5, .5, .5, 1});
    EXPECT_EQ(data->size(), n_particles - 2);
    EXPECT_EQ(data->getNDeactivated(), 2);
    EXPECT_EQ(*(data->end_positions() - 1), readdy::model::Vec3(.5, .5, .5));
}

TEST(ParticleData, markingForRemoval) {
    unsigned int n_particles = 100;
    auto &&data = std::make_unique<SingleCPUParticleData>(n_particles);
    EXPECT_EQ(0, data->size());
    for (auto &&i = 0; i < n_particles; i++) {
        data->addParticle({(double) i, (double) i, (double) i, 5});
    }
    EXPECT_EQ(n_particles, data->size());
    int n_deactivated = 0;
    for (size_t i = 0; i < 50; i++) {
        if (i % 3 == 0) {
            n_deactivated++;
            data->markForDeactivation(i);
        }
    }

    for (size_t i = 0; i < 50; i++) {
        if (i % 3 == 0) {
            EXPECT_TRUE(data->isMarkedForDeactivation(i));
        }
    }

    EXPECT_EQ(n_particles - n_deactivated, data->size());
    EXPECT_EQ(n_particles, data->getDeactivatedIndex());

    data->deactivateMarked();

    EXPECT_EQ(n_particles - n_deactivated, data->size());
}

TEST(ParticleData, empty) {
    auto data = std::make_unique<SingleCPUParticleData>(3);
    EXPECT_TRUE(data->empty());
    data->addParticle({4, 4, 4, 4});
    EXPECT_FALSE(data->empty());
}