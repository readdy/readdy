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
 * The observable Trajectory is a time series of all particles' positions. A single particle is represented by the
 * TrajectoryEntry. Unlike other observables, Trajectory must be constructed with a File object, because a trajectory
 * is not of much use except writing it to disk.
 *
 * @file Trajectory.h
 * @brief A trajectory keeps track of all particles' positions and saves them to a file.
 * @author chrisfroe
 * @author clonker
 * @date 12.12.16
 * @copyright BSD-3
 */

#pragma once

#include <array>
#include <memory>

#include <readdy/model/observables/Observable.h>

#include "TrajectoryEntry.h"

namespace readdy::model {
class Kernel;
namespace observables {

class Trajectory : public Observable<std::vector<TrajectoryEntry>> {
    using super = Observable<std::vector<TrajectoryEntry>>;
public:

    constexpr static auto &TRAJECTORY_GROUP_PATH = "/readdy/trajectory";

    Trajectory(model::Kernel *kernel, unsigned int stride);

    ~Trajectory() override;

    void evaluate() override;

    void flush() override;

    std::string_view type() const override;

protected:
    void initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

class FlatTrajectory : public Observable<std::vector<TrajectoryEntry>> {
    using super = Observable<std::vector<TrajectoryEntry>>;
public:

    FlatTrajectory(Kernel *kernel, unsigned int stride, bool useBlosc = true);

    ~FlatTrajectory() override;

    void evaluate() override;

    void flush() override;

    std::string_view type() const override;

protected:
    void initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;

    bool useBlosc{true};
};

}
}
