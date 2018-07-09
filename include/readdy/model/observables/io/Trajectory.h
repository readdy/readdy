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
 * The observable Trajectory is a time series of all particles' positions. A single particle is represented by the
 * TrajectoryEntry. Unlike other observables, Trajectory must be constructed with a File object, because a trajectory
 * is not of much use except writing it to disk.
 *
 * @file Trajectory.h
 * @brief A trajectory keeps track of all particles' positions and saves them to a file.
 * @author chrisfroe
 * @author clonker
 * @date 12.12.16
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <array>
#include <memory>

#include <readdy/model/observables/Observable.h>

#include "TrajectoryEntry.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
class Kernel;
NAMESPACE_BEGIN(observables)

class Trajectory : public Observable<std::vector<TrajectoryEntry>> {
    using super = Observable<std::vector<TrajectoryEntry>>;
public:

    const static std::string TRAJECTORY_GROUP_PATH;

    Trajectory(model::Kernel * kernel, unsigned int stride);

    ~Trajectory() override;

    void evaluate() override;

    void flush() override;

    std::string type() const override;

protected:
    void initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

class FlatTrajectory : public Observable<std::vector<TrajectoryEntry>> {
    using super = Observable<std::vector<TrajectoryEntry>>;
public:

    FlatTrajectory(Kernel* kernel, unsigned int stride);

    ~FlatTrajectory() override;

    FlatTrajectory(FlatTrajectory&&);

    void evaluate() override;

    void flush() override;

    std::string type() const override;

protected:
    void initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
