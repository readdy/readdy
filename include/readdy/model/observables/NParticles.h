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
 * @file NParticles.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include <vector>
#include "Observable.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class NParticles : public Observable<std::vector<unsigned long>> {

public:
    NParticles(Kernel* kernel, unsigned int stride);

    NParticles(Kernel* kernel, unsigned int stride, std::vector<std::string> typesToCount);

    NParticles(Kernel* kernel, unsigned int stride, std::vector<unsigned int> typesToCount);

    NParticles(const NParticles&) = delete;
    NParticles& operator=(const NParticles&) = delete;
    NParticles(NParticles&&) = default;
    NParticles& operator=(NParticles&&) = delete;

    void flush() override;

    virtual ~NParticles();

protected:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

    void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    std::vector<unsigned int> typesToCount;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)