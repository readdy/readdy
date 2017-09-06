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
 * @file PositionsObservable.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <vector>
#include <readdy/common/macros.h>
#include <readdy/model/Vec3.h>
#include "Observable.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class KernelContext;

class Positions : public Observable<std::vector<Vec3>> {
public:

    explicit Positions(Kernel* kernel, unsigned int stride = 1);

    void flush() override;

    Positions(Kernel* kernel, unsigned int stride, std::vector<std::string> typesToCount);

    Positions(Kernel* kernel, unsigned int stride, std::vector<unsigned int> typesToCount);

    Positions(const Positions&) = delete;
    Positions& operator=(const Positions&) = delete;
    Positions(Positions&&) = default;
    Positions& operator=(Positions&&) = delete;

    virtual ~Positions();

protected:

    void initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    std::vector<unsigned int> typesToCount;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
