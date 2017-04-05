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
 * @file TopologyReaction.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <vector>
#include <memory>
#include <functional>

#include <readdy/common/macros.h>

#include "Operation.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class GraphTopology;
NAMESPACE_BEGIN(reactions)

enum class Mode{
    RAISE, ROLLBACK
};

class TopologyReaction {
public:
    using reaction_operations = std::vector<op::OperationRef>;
    using rate_function = std::function<double(const GraphTopology &)>;

    TopologyReaction(reaction_operations operations, const rate_function &rate_function);

    double rate(const GraphTopology &topology) const;

    const reaction_operations &operations() const;

    const Mode& mode() const;

    Mode& mode();

private:
    rate_function rate_function_;
    reaction_operations operations_;
    Mode mode_ {Mode::RAISE};
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
