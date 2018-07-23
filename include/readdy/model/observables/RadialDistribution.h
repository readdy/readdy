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
 * @file RadialDistributionObservable.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright BSD-3
 */

#pragma once

#include <vector>
#include <readdy/common/macros.h>
#include <readdy/io/BloscFilter.h>
#include "Observable.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class RadialDistribution : public Observable<std::pair<std::vector<scalar>, std::vector<scalar>>> {
public:
    RadialDistribution(Kernel *kernel, stride_type stride, std::vector<scalar> binBorders,
                       std::vector<ParticleTypeId> typeCountFrom, std::vector<ParticleTypeId> typeCountTo,
                       scalar particleToDensity);

    RadialDistribution(Kernel *kernel, stride_type stride, const std::vector<scalar> &binBorders,
                       const std::vector<std::string> &typeCountFrom, const std::vector<std::string> &typeCountTo,
                       scalar particleToDensity);

    ~RadialDistribution() override;

    const std::vector<scalar> &getBinBorders() const;

    std::string type() const override;

    void evaluate() override;

    void flush() override;

protected:

    void setBinBorders(const std::vector<scalar> &binBorders);

    void initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
    std::vector<scalar> binBorders;
    std::vector<scalar> counts;
    std::vector<ParticleTypeId> typeCountFrom, typeCountTo;
    scalar particleToDensity;
    readdy::io::BloscFilter bloscFilter;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)