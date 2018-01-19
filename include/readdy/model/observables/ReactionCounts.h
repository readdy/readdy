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
#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/model/Particle.h>
#include <readdy/model/Context.h>
#include <readdy/model/reactions/Reaction.h>
#include <readdy/model/reactions/ReactionRecord.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class ReactionCounts : public Observable<reactions::reaction_counts_map> {
public:
    using reaction_counts_map = result_type;

    ReactionCounts(Kernel *kernel, stride_type stride);

    ReactionCounts(const ReactionCounts &) = delete;

    ReactionCounts &operator=(const ReactionCounts &) = delete;

    ReactionCounts(ReactionCounts &&) = default;

    ReactionCounts &operator=(ReactionCounts &&) = delete;

    virtual ~ReactionCounts();

    void flush() override;

    std::string type() const override;


protected:
    void initialize(Kernel *kernel) override;

    void initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) override;

    void append() override;

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
