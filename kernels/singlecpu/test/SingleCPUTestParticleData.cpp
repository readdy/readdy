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

struct ParticleDataTest : public ::testing::TestWithParam<bool> {

    std::unique_ptr<SingleCPUParticleData> createDataObject(std::size_t capacity) {
        return std::make_unique<SingleCPUParticleData>(capacity, GetParam());
    }
};

TEST_P(ParticleDataTest, initialization) {
    // everything should be deactivated
    auto data = createDataObject(10);

    for (size_t i = 0; i < 10; i++) {
        EXPECT_TRUE(data->isMarkedForDeactivation(i));
    }
}

TEST_P(ParticleDataTest, additionOfParticles) {
    // initial capacity of 3
    auto data = createDataObject(3);
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

TEST_P(ParticleDataTest, removalOfParticles) {
    unsigned int n_particles = 15;
    auto data = createDataObject(5);
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

TEST_P(ParticleDataTest, markingForRemoval) {
    unsigned int n_particles = 100;
    auto &&data = createDataObject(n_particles);
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

TEST_P(ParticleDataTest, markingForRemoval2) {
    unsigned int n_particles = 11;
    auto &&data = createDataObject(15);
    for (auto &&i = 0; i < n_particles; i++) {
        data->addParticle({(double) i, (double) i, (double) i, 5});
    }
    std::vector<boost::uuids::uuid> deactivatedIds = {
            *(data->begin_ids()+10),
            *(data->begin_ids()+7),
            *(data->begin_ids()+8),
            *(data->begin_ids()+0),
            *(data->begin_ids()+4)
    };
    data->markForDeactivation(10);
    data->markForDeactivation(7);
    data->markForDeactivation(8);
    data->markForDeactivation(0);
    data->markForDeactivation(4);
    data->deactivateMarked();
    EXPECT_EQ(data->size(), 6);

    bool noDeactivatedIds = true;
    for(auto it = deactivatedIds.begin(); it != deactivatedIds.end(); ++it) {
        noDeactivatedIds &= std::find(data->begin_ids(), data->end_ids(), *it) == data->end_ids();
    }
    EXPECT_TRUE(noDeactivatedIds);

    deactivatedIds = {
            *(data->begin_ids()+4),
            *(data->begin_ids()+2),
            *(data->begin_ids()+0),
            *(data->begin_ids()+3),
            *(data->begin_ids()+1)
    };
    data->markForDeactivation(4);
    data->markForDeactivation(2);
    data->markForDeactivation(0);
    data->markForDeactivation(3);
    data->markForDeactivation(1);
    data->deactivateMarked();
    noDeactivatedIds = true;
    for(auto it = deactivatedIds.begin(); it != deactivatedIds.end(); ++it) {
        noDeactivatedIds &= std::find(data->begin_ids(), data->end_ids(), *it) == data->end_ids();
    }
    EXPECT_TRUE(noDeactivatedIds);
    EXPECT_EQ(data->size(), 1);

    data->markForDeactivation(0);
    data->deactivateMarked();
    EXPECT_EQ(data->size(), 0);
}

TEST_P(ParticleDataTest, empty) {
    auto data = createDataObject(3);
    EXPECT_TRUE(data->empty());
    data->addParticle({4, 4, 4, 4});
    EXPECT_FALSE(data->empty());
}

INSTANTIATE_TEST_CASE_P(ParticleDataTests, ParticleDataTest, ::testing::Bool());