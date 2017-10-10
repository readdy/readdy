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
 * This file contains definitions of `getRecommendedTimeStep` and `getMaximumDisplacement` as well as
 * some utility with respect to validating type names.
 *
 * @file readdy/model/Utils.h
 * @brief Definitions of some utility for the readdy_model target
 * @author clonker
 * @date 17.11.16
 */

#pragma once
#include <readdy/model/Context.h>
#include <readdy/model/observables/Observable.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(util)

scalar getRecommendedTimeStep(unsigned int N, Context&);
scalar getMaximumDisplacement(Context&, scalar timeStep);

constexpr inline std::array<const char*, 8> invalidCharacterSequences() {
    return {{"[", "]", "(", ")", "->", "--", ":", "+"}};
};

static constexpr const char arrow[] = "->";
static constexpr const char bond[] = "--";

void validateTypeName(const std::string &typeName);

NAMESPACE_END(util)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
