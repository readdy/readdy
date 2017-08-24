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
 * @file TopologyFusionReaction.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include <readdy/common/common.h>
#include <readdy/common/ParticleTypeTuple.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(reactions)

class ExternalTopologyReaction {
public:
    ExternalTopologyReaction(const std::string &name, const util::particle_type_pair &types,
                             const util::particle_type_pair &types_to, scalar rate, scalar radius);

    ~ExternalTopologyReaction() = default;

    ExternalTopologyReaction(const ExternalTopologyReaction &) = default;

    ExternalTopologyReaction &operator=(const ExternalTopologyReaction &) = default;

    ExternalTopologyReaction(ExternalTopologyReaction &&) = default;

    ExternalTopologyReaction &operator=(ExternalTopologyReaction &&) = default;

    const std::string &name() const;

    const particle_type_type type1() const;

    const particle_type_type type2() const;

    const util::particle_type_pair &types() const;

    const particle_type_type type_to1() const;

    const particle_type_type type_to2() const;

    const util::particle_type_pair &types_to() const;

    const scalar rate() const;

    const scalar radius() const;

private:
    std::string _name;
    util::particle_type_pair _types;
    util::particle_type_pair _types_to;
    scalar _rate;
    scalar _radius;
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
