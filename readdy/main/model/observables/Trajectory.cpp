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
 * @file Trajectory.cpp
 * @brief Core library impl for Trajectory.h
 * @author chrisfroe
 * @author clonker
 * @date 12.12.16
 * @copyright GPL-3
 */

#include <readdy/model/Kernel.h>

#include <readdy/model/observables/io/Trajectory.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>
#include <readdy/model/observables/io/Types.h>

namespace readdy {
namespace model {
namespace observables {

const std::string Trajectory::TRAJECTORY_GROUP_PATH = "/readdy/trajectory";

struct Trajectory::Impl {
    std::unique_ptr<h5rd::VLENDataSet> dataSet {nullptr};
    std::unique_ptr<util::TimeSeriesWriter> time {nullptr};
    std::unique_ptr<util::CompoundH5Types> h5types {nullptr};
};


Trajectory::Trajectory(readdy::model::Kernel *const kernel, unsigned int stride)
        : super(kernel, stride), pimpl(std::make_unique<Impl>()) {
}

void Trajectory::evaluate() {
    result.clear();
    const auto &currentInput = kernel->stateModel().getParticles();
    std::for_each(currentInput.begin(), currentInput.end(), [this](const Particle &p) {
        result.emplace_back(p, kernel->context().particleTypes());
    });
}

void Trajectory::flush() {
    if (pimpl->dataSet) pimpl->dataSet->flush();
    if (pimpl->time) pimpl->time->flush();
}

Trajectory::~Trajectory() = default;

void Trajectory::initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) {
    pimpl->h5types = std::make_unique<util::CompoundH5Types>(util::getTrajectoryEntryTypes(file.parentFile()));
    h5rd::dimensions fs = {flushStride};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
    auto group = file.createGroup(
            std::string(TRAJECTORY_GROUP_PATH + (dataSetName.length() > 0 ? "/" + dataSetName : "")));
    pimpl->dataSet = group.createVLENDataSet("records", fs, dims,
                                             std::get<0>(*pimpl->h5types), std::get<1>(*pimpl->h5types));
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void Trajectory::append() {
    pimpl->dataSet->append({1}, &result);
    pimpl->time->append(t_current);
}

std::string Trajectory::type() const {
    return "Trajectory";
}

struct FlatTrajectory::Impl {
    std::unique_ptr<h5rd::DataSet> dataSet {nullptr};
    std::unique_ptr<h5rd::DataSet> limits {nullptr};
    std::unique_ptr<util::TimeSeriesWriter> time {nullptr};
    std::unique_ptr<util::CompoundH5Types> h5types {nullptr};
    std::size_t current_limits[2]{0, 0};
};

FlatTrajectory::FlatTrajectory(Kernel *const kernel, unsigned int stride)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {}

void FlatTrajectory::initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->dataSet) {
        pimpl->h5types = std::make_unique<util::CompoundH5Types>(util::getTrajectoryEntryTypes(file.parentFile()));
        auto group = file.createGroup(
                std::string(Trajectory::TRAJECTORY_GROUP_PATH + (dataSetName.length() > 0 ? "/" + dataSetName : "")));
        {
            h5rd::dimensions fs = {flushStride};
            h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
            io::BloscFilter filter;
            pimpl->dataSet = group.createDataSet("records", fs, dims, std::get<0>(*pimpl->h5types),
                                                 std::get<1>(*pimpl->h5types), {&filter});
        }
        {
            h5rd::dimensions fs = {flushStride, 2};
            h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS, 2};
            pimpl->limits = group.createDataSet<std::size_t>("limits", fs, dims);
        }
        pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);

    }
}

void FlatTrajectory::evaluate() {
    result.clear();
    const auto &currentInput = kernel->stateModel().getParticles();
    std::for_each(currentInput.begin(), currentInput.end(), [this](const Particle &p) {
        result.emplace_back(p, kernel->context().particleTypes());
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

std::string FlatTrajectory::type() const {
    return "FlatTrajectory";
}

FlatTrajectory::FlatTrajectory(FlatTrajectory &&) = default;

FlatTrajectory::~FlatTrajectory() = default;
}
}
}