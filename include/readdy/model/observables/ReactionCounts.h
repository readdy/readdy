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
 * @file ReactionCounts.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "Observable.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class ReactionCounts : public Observable<std::tuple<std::vector<std::size_t>, std::vector<std::size_t>>> {
public:
    ReactionCounts(Kernel *const kernel, unsigned int stride);

    virtual ~ReactionCounts();

    virtual void flush() override;

protected:
    virtual void initialize(Kernel *const kernel) override;

    virtual void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    virtual void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};


NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
