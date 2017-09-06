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
 * @file TimeSeriesWriter.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once


#include <h5rd/h5rd.h>

#include <readdy/common/common.h>
#include <readdy/io/BloscFilter.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)
NAMESPACE_BEGIN(util)

class TimeSeriesWriter {
public:
    TimeSeriesWriter(h5rd::Group &group, unsigned int chunkSize, const std::string &dsName = "time")
            : dataSet(group.createDataSet<time_step_type>(dsName, {chunkSize}, {h5rd::UNLIMITED_DIMS}, {&bloscFilter})) {}

    ~TimeSeriesWriter() = default;

    TimeSeriesWriter(const TimeSeriesWriter &) = delete;

    TimeSeriesWriter &operator=(const TimeSeriesWriter &) = delete;

    TimeSeriesWriter(TimeSeriesWriter &&) = default;

    TimeSeriesWriter &operator=(TimeSeriesWriter &&) = default;

    void append(const time_step_type t) {
        dataSet->append({1}, &t);
    }

    void flush() {
        dataSet->flush();
    }

private:
    io::BloscFilter bloscFilter {};
    std::unique_ptr<h5rd::DataSet> dataSet;
};

NAMESPACE_END(util)
NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)