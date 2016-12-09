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


/**
 * << detailed description >>
 *
 * @file TestIO.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12.12.16
 * @copyright GNU Lesser General Public License v3.0
 */

#include <numeric>

#include <gtest/gtest.h>

#include <readdy/io/File.h>

namespace {

using file = readdy::io::File;

struct TestIO : public ::testing::Test {
protected:
    std::string tempFile;

    void SetUp() override {
        char out[4096];
#if READDY_WINDOWS
    if (GetTempFileName(".", "rdy", 0, out) != 0) {
		// Strip the ".\\" prefix
		if (strncmp(out, ".\\", 2) == 0)
		{
			memmove(out, out + 2, strlen(out));
		}
		// Create file
		fclose(fopen(out, "w"));
	}
#else
        close(mkstemp(out));
#endif
        tempFile = std::string(out);
    }

    void TearDown() override {
        std::remove(tempFile.c_str());
    }
};

TEST_F(TestIO, TestGroupWrite) {
    std::vector<double> data;
    data.resize(100);
    std::iota(data.begin(), data.end(), 0);
    {
        file f(tempFile, file::Action::CREATE, file::Flag::OVERWRITE);
        f.write("foo", data);
    }
}

}