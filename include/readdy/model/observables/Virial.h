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


/**
 * << detailed description >>
 *
 * @file Virial.h
 * @brief << brief description >>
 * @author clonker
 * @date 1/17/18
 */


#pragma once

#include <readdy/common/macros.h>
#include "Observable.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class Virial : public Observable<Matrix33> {
    using super = Observable<Matrix33>;
public:

    Virial(Kernel* kernel, stride_type stride);

    Virial(Virial &&) = default;
    Virial& operator=(Virial &&) = default;
    Virial(const Virial &) = delete;
    Virial &operator=(const Virial &) = delete;

    ~Virial();

    std::string type() const override;

    void flush() override;

protected:
    void initialize(Kernel * kernel) override;

    void initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
