/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file TimeSeriesWriter.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright BSD-3
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