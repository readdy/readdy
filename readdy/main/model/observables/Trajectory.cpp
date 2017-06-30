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
 * @file Trajectory.cpp
 * @brief Core library impl for Trajectory.h
 * @author chrisfroe
 * @author clonker
 * @date 12.12.16
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/io/Trajectory.h>
#include <readdy/io/File.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>
#include <readdy/model/observables/io/Types.h>

namespace io = readdy::io;

using particle_t = readdy::model::Particle;
using vec_t = readdy::model::Vec3;
namespace readdy {
namespace model {
namespace observables {

const std::string Trajectory::TRAJECTORY_GROUP_PATH = "/readdy/trajectory";

struct Trajectory::Impl {
    std::unique_ptr<io::VLENDataSet> dataSet;
    std::unique_ptr<util::TimeSeriesWriter> time;
};


Trajectory::Trajectory(readdy::model::Kernel *const kernel, unsigned int stride)
        : super(kernel, stride), pimpl(std::make_unique<Impl>()) {
}

void Trajectory::evaluate() {
    result.clear();
    const auto &currentInput = kernel->getKernelStateModel().getParticles();
    std::for_each(currentInput.begin(), currentInput.end(), [this](const Particle &p) {
        result.push_back(TrajectoryEntry{p});
    });
}

void Trajectory::flush() {
    if (pimpl->dataSet) pimpl->dataSet->flush();
    if (pimpl->time) pimpl->time->flush();

}

Trajectory::~Trajectory() = default;

void Trajectory::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->dataSet) {
        std::vector<readdy::io::h5::dims_t> fs = {flushStride};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        auto group = file.createGroup(
                std::string(TRAJECTORY_GROUP_PATH + (dataSetName.length() > 0 ? "/" + dataSetName : "")));
        auto dataSet = std::make_unique<io::VLENDataSet>(group.createVLENDataSet(
                "records", fs, dims, util::TrajectoryEntryMemoryType(), util::TrajectoryEntryFileType()));
        pimpl->dataSet = std::move(dataSet);
        pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
    }
}

void Trajectory::append() {
    pimpl->dataSet->append({1}, &result);
    pimpl->time->append(t_current);
}

struct FlatTrajectory::Impl {
    std::unique_ptr<readdy::io::DataSet> dataSet;
    std::unique_ptr<readdy::io::DataSet> limits;
    std::unique_ptr<util::TimeSeriesWriter> time;
    std::size_t current_limits[2]{0, 0};
};

FlatTrajectory::FlatTrajectory(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride),
                                                                            pimpl(std::make_unique<Impl>()) {}

void FlatTrajectory::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->dataSet) {
        auto group = file.createGroup(
                std::string(Trajectory::TRAJECTORY_GROUP_PATH + (dataSetName.length() > 0 ? "/" + dataSetName : "")));
        {
            std::vector<readdy::io::h5::dims_t> fs = {flushStride};
            std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
            auto dataSet = std::make_unique<readdy::io::DataSet>(group.createDataSet(
                    "records", fs, dims, util::TrajectoryEntryMemoryType(), util::TrajectoryEntryFileType()
            ));
            pimpl->dataSet = std::move(dataSet);
        }
        {
            std::vector<readdy::io::h5::dims_t> fs = {flushStride, 2};
            std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, 2};
            auto limits = std::make_unique<readdy::io::DataSet>(group.createDataSet<std::size_t>("limits", fs, dims));
            pimpl->limits = std::move(limits);
        }
        pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);

    }
}

void FlatTrajectory::evaluate() {
    result.clear();
    const auto &currentInput = kernel->getKernelStateModel().getParticles();
    std::for_each(currentInput.begin(), currentInput.end(), [this](const Particle &p) {
        result.push_back(TrajectoryEntry{p});
    });
}

void FlatTrajectory::flush() {
    if (pimpl->dataSet) pimpl->dataSet->flush();
    if (pimpl->time) pimpl->time->flush();
    if (pimpl->limits) pimpl->limits->flush();
}

void FlatTrajectory::append() {
    pimpl->current_limits[0] = pimpl->current_limits[1];
    pimpl->current_limits[1] += result.size();
    pimpl->dataSet->append({result.size()}, result.data());
    pimpl->time->append(t_current);
    pimpl->limits->append({1, 2}, pimpl->current_limits);
}

FlatTrajectory::FlatTrajectory(FlatTrajectory &&) = default;

FlatTrajectory::~FlatTrajectory() = default;
}
}
}