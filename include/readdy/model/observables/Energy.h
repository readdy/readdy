/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * @file Energy.h
 * @brief Definitions of the energy observable
 * @author clonker
 * @date 12/21/17
 */


#pragma once


#include <readdy/common/common.h>
#include "Observable.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class Energy : public Observable<scalar> {
public:
    Energy(Kernel *kernel, stride_type stride);

    Energy(const Energy &) = delete;

    Energy &operator=(const Energy &) = delete;

    Energy(Energy &&) = delete;

    Energy &operator=(Energy &&) = delete;

    ~Energy() override;

    void flush() override;

    void evaluate() override;

    std::string type() const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

    void initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) override;

    void append() override;

};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
