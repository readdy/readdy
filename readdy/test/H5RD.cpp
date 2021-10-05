/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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

/**
 * << detailed description >>
 *
 * @file TestMain.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 04.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <utility>
#include <iostream>

#include <h5rd/h5rd.h>
#include <catch2/catch.hpp>

TEST_CASE("Sanity checks", "[sanity]") {
    using namespace h5rd;

    {
        auto f = File::create("test.h5", File::Flag::OVERWRITE);
        std::cout << "/ exists: " << f->exists("/") << std::endl;
        auto g = f->createGroup("/foo/bar");
        std::cout << "foo/bar exists: " << f->exists("/foo/bar") << std::endl;

        g.write("miau", "blubs");
        std::vector<int> zahlen{1, 2, 3, 4, 5, 6, 7};
        g.write("zahlen", zahlen);
        std::vector<int> readZahlen;
        g.read("zahlen", readZahlen);
        for (auto x : readZahlen) {
            std::cout << x;
        }
        std::cout << std::endl;

        f->close();

        auto f2 = File::create("test.h5", File::Flag::OVERWRITE);
        std::cout << "/ exists: " << f2->exists("/") << std::endl;
        auto g2 = f2->createGroup("/foo/bar");
        std::cout << "foo/bar exists: " << f2->exists("/foo/bar") << std::endl;
    }
}

struct Stuff {
    int a;
    float x;
    std::array<short, 3> xyz;

    bool operator==(const Stuff &rhs) const {
        return a == rhs.a && x == rhs.x && xyz == rhs.xyz;
    }
};

auto getCompoundTypes(h5rd::Object::ParentFileRef parentFile) {
    using namespace h5rd;
    NativeCompoundType stuffType = NativeCompoundTypeBuilder(sizeof(Stuff), std::move(parentFile))
            .insert<decltype(std::declval<Stuff>().a)>("a", offsetof(Stuff, a))
            .insert<decltype(std::declval<Stuff>().x)>("x", offsetof(Stuff, x))
            .insertStdArray<decltype(std::declval<Stuff>().xyz)>("xyz", offsetof(Stuff, xyz))
            .build();
    STDCompoundType stdStuffType(stuffType);
    return std::make_tuple(std::move(stuffType), std::move(stdStuffType));
}

TEST_CASE("Check reading and writing of compound types", "[compound]") {
    using namespace h5rd;
    {
        auto f = File::create("test.h5", File::Flag::OVERWRITE);
        auto group = f->createGroup("/my/compound/group");

        Stuff s1{3, 5.f, {{1, 1, 1}}};
        Stuff s2{4, 5.f, {{2, 2, 2}}};
        Stuff s3{6, 7.f, {{3, 3, 3}}};
        std::vector<Stuff> stuffs{s1, s2, s3};

        auto types = getCompoundTypes(f->parentFile());

        auto ds = group.createDataSet("stuffs", {3}, {UNLIMITED_DIMS}, std::get<0>(types), std::get<1>(types));
        ds->append({stuffs.size()}, stuffs.data());

    }
    {
        auto f = File::open("test.h5", File::Flag::READ_ONLY);
        std::cout << f->subgroups()[0] << '\n';
        auto group = f->getSubgroup("/my/compound/group");
        auto types = getCompoundTypes(f->parentFile());
        std::vector<Stuff> stuffs;
        group.read("stuffs", stuffs, &std::get<0>(types), &std::get<1>(types));

        REQUIRE(stuffs.size() == 3);
        for (const auto &s : stuffs) {
            std::cout << "stuff: a=" << s.a << ", x=" << s.x << ", xyz=" << s.xyz.at(0) << "," << s.xyz.at(1) << ","
                      << s.xyz.at(2) << std::endl;
        }

        std::vector<Stuff> stridedStuff;
        group.read("stuffs", stridedStuff, &std::get<0>(types), &std::get<1>(types), {2});

        REQUIRE(stridedStuff.size() == 2);

        REQUIRE(stridedStuff.at(0).a == 3);
        REQUIRE(stridedStuff.at(1).a == 6);

        std::vector<Stuff> selectedStuff;
        group.readSelection("stuffs", selectedStuff, &std::get<0>(types), &std::get<1>(types), {1}, {1}, {2});
        REQUIRE(selectedStuff.size() == 2);

        REQUIRE(selectedStuff.at(0).a == 4);
        REQUIRE(selectedStuff.at(1).a == 6);
    }

}

TEST_CASE("Check reading and writing of compound types with repeated opening and closing.", "[compound]") {
    using namespace h5rd;
    auto f = File::create("test.h5", File::Flag::OVERWRITE);
    auto group = f->createGroup("/my/compound/group");

    std::vector<Stuff> stuffs;
    {
        Stuff s1{1, 1.f, {{1, 1, 1}}};
        Stuff s2{2, 2.f, {{2, 2, 2}}};
        Stuff s3{3, 3.f, {{3, 3, 3}}};
        stuffs.push_back(s1);
        stuffs.push_back(s2);
        stuffs.push_back(s3);
    }

    auto types = getCompoundTypes(f->parentFile());


    auto ds = group.createDataSet("stuffs", {3}, {UNLIMITED_DIMS}, std::get<0>(types), std::get<1>(types));
    {
        ds->append({stuffs.size()}, stuffs.data());
        ds->append({stuffs.size()}, stuffs.data());
        f->close();
        REQUIRE(f->closed());

        f.reset();
    }

    {
        auto f2 = File::open("test.h5", File::Flag::READ_ONLY);
        auto group2 = f2->getSubgroup("/my/compound/group");
        auto types2 = getCompoundTypes(f2->parentFile());
        std::vector<Stuff> stuffs2;
        group2.read("stuffs", stuffs2, &std::get<0>(types2), &std::get<1>(types2));

        REQUIRE(stuffs2.size() == 6);

        int i = 1;
        for (const auto &s : stuffs2) {
            REQUIRE(s.a == 1 + ((i - 1) % 3));
            REQUIRE(s.x == 1 + ((i - 1) % 3));
            REQUIRE(s.xyz.at(0) == 1 + ((i - 1) % 3));
            REQUIRE(s.xyz.at(1) == 1 + ((i - 1) % 3));
            REQUIRE(s.xyz.at(2) == 1 + ((i - 1) % 3));
            ++i;
        }
    }

    auto f3 = File::create("test.h5", File::Flag::OVERWRITE);
    auto gtest123 = f3->createGroup("/test123");

    gtest123.write("hier", "stehtwas");

    auto types123 = getCompoundTypes(f3->parentFile());
    auto ds213 = gtest123.createDataSet("stuffs", {3}, {UNLIMITED_DIMS}, std::get<0>(types123), std::get<1>(types123));
    ds.reset();
    ds213->append({stuffs.size()}, stuffs.data());


}

TEST_CASE("Reading and writing of VLEN data set", "[vlen]") {
    using namespace h5rd;
    auto f = File::create("test.h5", File::Flag::OVERWRITE);
    auto group = f->createGroup("/my/vlen/group");

    auto ds = group.createVLENDataSet<double>("myds", {3}, {UNLIMITED_DIMS});

    {
        std::vector<std::vector<double>> data{{5., 7.},
                                              {0., 1., 2., 3.}};
        ds->append(data);
    }

    {
        std::vector<std::vector<double>> data;
        group.readVLEN("myds", data);

        REQUIRE(data.size() == 2);
        REQUIRE(data.at(0).size() == 2);
        REQUIRE(data.at(0).at(0) == 5);
        REQUIRE(data.at(0).at(1) == 7);
        REQUIRE(data.at(1).size() == 4);
        REQUIRE(data.at(1).at(0) == 0);
        REQUIRE(data.at(1).at(1) == 1);
        REQUIRE(data.at(1).at(2) == 2);
        REQUIRE(data.at(1).at(3) == 3);
    }
}

TEST_CASE("Reading and writing of compound VLEN data set", "[vlen]") {
    using namespace h5rd;
    auto f = File::create("test.h5", File::Flag::OVERWRITE);
    auto group = f->createGroup("/my/vlen/compound/group");

    std::vector<std::vector<Stuff>> allStuffs(3);
    {
        auto &stuffs = allStuffs.at(0);
        Stuff s1{1, 1.f, {{1, 1, 1}}};
        Stuff s2{2, 2.f, {{2, 2, 2}}};
        Stuff s3{3, 3.f, {{3, 3, 3}}};
        stuffs.push_back(s1);
        stuffs.push_back(s2);
        stuffs.push_back(s3);
    }
    {
        auto &stuffs = allStuffs.at(1);
        Stuff s2{2, 2.f, {{5, 5, 5}}};
        Stuff s3{3, 3.f, {{6, 6, 6}}};
        stuffs.push_back(s2);
        stuffs.push_back(s3);
    }
    {
        auto &stuffs = allStuffs.at(2);
        Stuff s1{1, 1.f, {{1, 1, 1}}};
        Stuff s2{2, 2.f, {{5, 5, 5}}};
        Stuff s3{3, 3.f, {{6, 6, 6}}};
        Stuff s4{4, 4.f, {{7, 7, 7}}};
        stuffs.push_back(s1);
        stuffs.push_back(s2);
        stuffs.push_back(s3);
        stuffs.push_back(s4);
    }

    auto types = getCompoundTypes(f->parentFile());


    auto ds = group.createVLENDataSet("stuffs", {3}, {UNLIMITED_DIMS}, std::get<0>(types), std::get<1>(types));
    {
        ds->append(allStuffs);
    }

    {
        {
            auto fs = group.getDataset("stuffs", std::get<0>(types), std::get<1>(types))->getFileSpace();
            REQUIRE(fs->dims().size() == 1);
            REQUIRE(fs->dims()[0] == 3);
        }
        std::vector<std::vector<Stuff>> readBackStuff;
        group.readVLEN("stuffs", readBackStuff, &std::get<0>(types), &std::get<1>(types));
        REQUIRE(3 == readBackStuff.size());
        REQUIRE(3 == readBackStuff.at(0).size());
        for (int i = 0; i < 3; ++i) {
            REQUIRE(readBackStuff.at(0).at(i).x == i + 1);
            REQUIRE(readBackStuff.at(0).at(i).a == i + 1);
            for (int j = 0; j < 3; ++j) {
                REQUIRE(readBackStuff.at(0).at(i).xyz.at(j) == i + 1);
            }
        }
        REQUIRE(2 == readBackStuff.at(1).size());
        for (int i = 0; i < 2; ++i) {
            REQUIRE(readBackStuff.at(1).at(i).x == i + 2);
            REQUIRE(readBackStuff.at(1).at(i).a == i + 2);
            for (int j = 0; j < 3; ++j) {
                REQUIRE(readBackStuff.at(1).at(i).xyz.at(j) == i + 5);
            }
        }
    }
    {
        std::vector<std::vector<Stuff>> readBackStuff;
        group.readVLENSelection("stuffs", readBackStuff, &std::get<0>(types), &std::get<1>(types), {0}, {1}, {1});
        REQUIRE(1 == readBackStuff.size());
        REQUIRE(allStuffs[0].size() == readBackStuff[0].size());
        REQUIRE(std::equal(readBackStuff[0].begin(), readBackStuff[0].end(), allStuffs[0].begin()));
    }
    {
        std::vector<std::vector<Stuff>> readBackStuff;
        group.readVLENSelection("stuffs", readBackStuff, &std::get<0>(types), &std::get<1>(types), {1}, {1}, {1});
        REQUIRE(1 == readBackStuff.size());
        REQUIRE(allStuffs[1].size() == readBackStuff[0].size());
        REQUIRE(std::equal(readBackStuff[0].begin(), readBackStuff[0].end(), allStuffs[1].begin()));
    }
    {
        std::vector<std::vector<Stuff>> readBackStuff;
        group.readVLENSelection("stuffs", readBackStuff, &std::get<0>(types), &std::get<1>(types), {1}, {1}, {2});
        REQUIRE(2 == readBackStuff.size());
        REQUIRE(allStuffs[1].size() == readBackStuff[0].size());
        REQUIRE(allStuffs[2].size() == readBackStuff[1].size());
        REQUIRE(std::equal(readBackStuff[0].begin(), readBackStuff[0].end(), allStuffs[1].begin()));
        REQUIRE(std::equal(readBackStuff[1].begin(), readBackStuff[1].end(), allStuffs[2].begin()));
    }

    {
        REQUIRE(1 == f->subgroups().size());
        REQUIRE(1 == group.containedDataSets().size());
        REQUIRE("stuffs" == group.containedDataSets()[0]);
    }
}
