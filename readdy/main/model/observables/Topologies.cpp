/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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


#include <readdy/model/observables/Topologies.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>
#include <readdy/model/observables/io/Types.h>

namespace readdy {
namespace model {
namespace observables {

struct Topologies::Impl {
    std::array<std::size_t, 2> currentLimitsParticles {{0_z, 0_z}};
    std::unique_ptr<h5rd::DataSet> dataSetParticles {nullptr};
    std::unique_ptr<h5rd::DataSet> limitsParticles {nullptr};

    std::array<std::size_t, 2> currentLimitsEdges {{0_z, 0_z}};
    std::unique_ptr<h5rd::DataSet> dataSetEdges {nullptr};
    std::unique_ptr<h5rd::DataSet> limitsEdges {nullptr};

    std::unique_ptr<util::TimeSeriesWriter> time {nullptr};
};

Topologies::Topologies(Kernel *kernel, stride_type stride)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {}

void Topologies::evaluate() {
    result.clear();
    for(auto topologyPtr : kernel->stateModel().getTopologies()) {
        top::TopologyRecord record;

        record.particleIndices = topologyPtr->getParticles();
        kernel->stateModel().toDenseParticleIndices(record.particleIndices.begin(), record.particleIndices.end());

        for(auto&& edge : topologyPtr->graph().edges()) {
            record.edges.emplace_back(std::get<0>(edge)->particleIndex, std::get<1>(edge)->particleIndex);
        }

        result.push_back(record);
    }
}

void Topologies::flush() {
    if(pimpl->dataSetParticles) pimpl->dataSetParticles->flush();
    if(pimpl->limitsParticles) pimpl->limitsParticles->flush();
    if(pimpl->dataSetEdges) pimpl->dataSetEdges->flush();
    if(pimpl->limitsEdges) pimpl->limitsEdges->flush();
    if(pimpl->time) pimpl->time->flush();
}

std::string Topologies::type() const {
    return "Topologies";
}

void Topologies::initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) {
    auto group = file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName);
    io::BloscFilter filter;
    {
        // data
        h5rd::dimensions fs = {flushStride};
        h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
        pimpl->dataSetParticles = group.createDataSet<std::size_t>("particles", fs, dims, {&filter});
        fs = {flushStride, 2};
        dims = {h5rd::UNLIMITED_DIMS, 2};
        pimpl->dataSetEdges = group.createDataSet<std::size_t>("edges", fs, dims, {&filter});
    }
    {
        // limits
        h5rd::dimensions fs = {flushStride, 2};
        h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS, 2};
        pimpl->limitsParticles = group.createDataSet<std::size_t>("limitsParticles", fs, dims);
        pimpl->limitsEdges = group.createDataSet<std::size_t>("limitsEdges", fs, dims);
    }

    pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void Topologies::append() {
    std::size_t totalNParticles {0}, totalNEdges {0};
    for(const auto &record : result) {
        totalNParticles += record.particleIndices.size();
        totalNEdges += record.edges.size();
    }
    // advance limits by total number of particles in all topologies + #topologies for the prefix
    pimpl->currentLimitsParticles[0] = pimpl->currentLimitsParticles[1];
    pimpl->currentLimitsParticles[1] += totalNParticles + result.size();

    // advance limits by total number of edges in all topologies + #topologies for the prefix
    pimpl->currentLimitsEdges[0] = pimpl->currentLimitsEdges[1];
    pimpl->currentLimitsEdges[1] += totalNEdges + result.size();

    std::vector<std::size_t> flatParticles;
    flatParticles.reserve(totalNParticles + result.size());
    std::vector<std::array<std::size_t, 2>> flatEdges;
    flatEdges.reserve(totalNEdges + result.size());

    for(const auto &r : result) {
        flatParticles.push_back(r.particleIndices.size());
        flatParticles.insert(std::end(flatParticles), std::begin(r.particleIndices), std::end(r.particleIndices));
        flatEdges.push_back(std::array<std::size_t, 2>{{r.edges.size(), 0}});
        for(const auto &edge : r.edges) {
            flatEdges.push_back(std::array<std::size_t, 2>{{std::get<0>(edge), std::get<1>(edge)}});
        }
    }

    pimpl->dataSetParticles->append({flatParticles.size()}, flatParticles.data());
    pimpl->limitsParticles->append({1, 2}, pimpl->currentLimitsParticles.data());
    pimpl->dataSetEdges->append({flatEdges.size(), 2}, &flatEdges[0][0]);
    pimpl->limitsEdges->append({1, 2}, pimpl->currentLimitsEdges.data());
    pimpl->time->append(t_current);
}

Topologies::~Topologies() = default;

}
}
}
