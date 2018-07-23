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
 * @file ParticlesObservable.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GPL-3
 */

#include <readdy/model/observables/Particles.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy {
namespace model {
namespace observables {


struct Particles::Impl {
    std::unique_ptr<h5rd::VLENDataSet> dataSetTypes;
    std::unique_ptr<h5rd::VLENDataSet> dataSetIds;
    std::unique_ptr<h5rd::VLENDataSet> dataSetPositions;
    std::unique_ptr<util::TimeSeriesWriter> time;
};

Particles::Particles(Kernel *const kernel, stride_type stride) : Observable(kernel, stride),
                                                                  pimpl(std::make_unique<Impl>()) {}

void Particles::initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) {
    using particle_t = readdy::model::Particle;
    h5rd::dimensions fs = {flushStride};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
    auto group = file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName);
    pimpl->dataSetTypes = group.createVLENDataSet<ParticleTypeId>("types", fs, dims);
    pimpl->dataSetIds = group.createVLENDataSet<particle_t::id_type>("ids", fs, dims);
    pimpl->dataSetPositions = group.createVLENDataSet("positions", fs, dims,
                                                      h5rd::NativeArrayDataSetType<scalar, 3>(group.parentFile()),
                                                      h5rd::STDArrayDataSetType<scalar, 3>(group.parentFile()));
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void Particles::append() {
    {
        auto &types = std::get<0>(result);
        log::debug("appending {} types ", types.size());
        for(auto t : types) {
            log::debug("    -> {}", t);
        }
        pimpl->dataSetTypes->append({1}, &types);
    }
    {
        auto &ids = std::get<1>(result);
        pimpl->dataSetIds->append({1}, &ids);
    }
    {
        pimpl->dataSetPositions->append({1}, &std::get<2>(result));
    }
    pimpl->time->append(t_current);
}

void Particles::flush() {
    if (pimpl->dataSetTypes) pimpl->dataSetTypes->flush();
    if (pimpl->dataSetIds) pimpl->dataSetIds->flush();
    if (pimpl->dataSetPositions) pimpl->dataSetPositions->flush();
    if (pimpl->time) pimpl->time->flush();
}

std::string Particles::type() const {
    return "Particles";
}

Particles::~Particles() = default;


}
}
}
