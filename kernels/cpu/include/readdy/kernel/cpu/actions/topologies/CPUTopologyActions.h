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
 * @file CPUTopologyActions.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.02.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include <readdy/model/topologies/Topology.h>
#include <readdy/model/topologies/potentials/TopologyPotentialActions.h>
#include <readdy/kernel/cpu/CPUStateModel.h>
#include <readdy/model/topologies/GraphTopology.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(kernel)
NAMESPACE_BEGIN(cpu)
NAMESPACE_BEGIN(actions)
NAMESPACE_BEGIN(top)

class CPUCalculateHarmonicBondPotential : public readdy::model::top::pot::CalculateHarmonicBondPotential {

    const harmonic_bond *const potential;
    CPUStateModel::data_type *const data;

public:
    CPUCalculateHarmonicBondPotential(const readdy::model::Context *context,
                                      CPUStateModel::data_type *data,
                                      const harmonic_bond *potential);

    readdy::scalar perform(const readdy::model::top::Topology *topology) override;
};


class CPUCalculateHarmonicAnglePotential : public readdy::model::top::pot::CalculateHarmonicAnglePotential {
    const harmonic_angle *const potential;
    CPUStateModel::data_type *const data;
public:
    CPUCalculateHarmonicAnglePotential(const readdy::model::Context *context, CPUStateModel::data_type *data,
                                       const harmonic_angle *potential);

    readdy::scalar perform(const readdy::model::top::Topology *topology) override;
};

class CPUCalculateCosineDihedralPotential : public readdy::model::top::pot::CalculateCosineDihedralPotential {
    const cos_dihedral *const potential;
    CPUStateModel::data_type *const data;
public:
    CPUCalculateCosineDihedralPotential(const readdy::model::Context *context,
                                        CPUStateModel::data_type *data,
                                        const cos_dihedral *pot);

    readdy::scalar perform(const readdy::model::top::Topology *topology) override;
};

NAMESPACE_BEGIN(reactions)
NAMESPACE_BEGIN(op)

class CPUChangeParticleType : public readdy::model::top::reactions::actions::ChangeParticleType {
    CPUStateModel::data_type *const data;
public:
    CPUChangeParticleType(CPUStateModel::data_type *data, model::top::GraphTopology *topology, const vertex &v,
                          const particle_type_type &type_to);

    void execute() override;

    void undo() override;

};

NAMESPACE_END(op)
NAMESPACE_END(reactions)

NAMESPACE_END(top)
NAMESPACE_END(actions)
NAMESPACE_END(cpu)
NAMESPACE_END(kernel)
NAMESPACE_END(readdy)
