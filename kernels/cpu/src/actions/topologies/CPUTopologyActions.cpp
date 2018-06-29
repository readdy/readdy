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
 * << detailed description >>
 *
 * @file CPUTopologyActions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12/21/17
 */

#include "readdy/kernel/cpu/actions/topologies/CPUTopologyActions.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace top {


CPUCalculateHarmonicBondPotential::CPUCalculateHarmonicBondPotential(const model::Context *const context,
                                                                     CPUStateModel::data_type *const data,
                                                                     const harmonic_bond *const potential)
        : CalculateHarmonicBondPotential(context), potential(potential), data(data) {}

readdy::scalar CPUCalculateHarmonicBondPotential::perform(const readdy::model::top::Topology *const topology) {
    scalar energy = 0;
    const auto &particleIndices = topology->getParticles();
    for (const auto &bond : potential->getBonds()) {
        Vec3 forceUpdate{0, 0, 0};
        auto &e1 = data->entry_at(particleIndices.at(bond.idx1));
        auto &e2 = data->entry_at(particleIndices.at(bond.idx2));
        const auto x_ij = bcs::shortestDifference(e1.pos, e2.pos, context->boxSize(),
                                                  context->periodicBoundaryConditions());
        potential->calculateForce(forceUpdate, x_ij, bond);
        e1.force += forceUpdate;
        e2.force -= forceUpdate;
        energy += potential->calculateEnergy(x_ij, bond);
    }
    return energy;
}


CPUCalculateHarmonicAnglePotential::CPUCalculateHarmonicAnglePotential(const model::Context *const context,
                                                                       CPUStateModel::data_type *const data,
                                                                       const harmonic_angle *const potential)
        : CalculateHarmonicAnglePotential(context), potential(potential), data(data) {}

readdy::scalar CPUCalculateHarmonicAnglePotential::perform(const readdy::model::top::Topology *const topology) {
    scalar energy = 0;
    const auto &particleIndices = topology->getParticles();


    for (const auto &angle : potential->getAngles()) {
        auto &e1 = data->entry_at(particleIndices.at(angle.idx1));
        auto &e2 = data->entry_at(particleIndices.at(angle.idx2));
        auto &e3 = data->entry_at(particleIndices.at(angle.idx3));
        const auto x_ji = bcs::shortestDifference(e2.pos, e1.pos, context->boxSize(),
                                                  context->periodicBoundaryConditions());
        const auto x_jk = bcs::shortestDifference(e2.pos, e3.pos, context->boxSize(),
                                                  context->periodicBoundaryConditions());
        energy += potential->calculateEnergy(x_ji, x_jk, angle);
        potential->calculateForce(e1.force, e2.force, e3.force, x_ji, x_jk, angle);
    }
    return energy;
}


CPUCalculateCosineDihedralPotential::CPUCalculateCosineDihedralPotential(const model::Context *const context,
                                                                         CPUStateModel::data_type *const data,
                                                                         const cos_dihedral *const pot)
        : CalculateCosineDihedralPotential(context), potential(pot), data(data) {
}

readdy::scalar CPUCalculateCosineDihedralPotential::perform(const readdy::model::top::Topology *const topology) {
    scalar energy = 0;
    const auto &particleIndices = topology->getParticles();

    for (const auto &dih : potential->getDihedrals()) {
        auto &e_i = data->entry_at(particleIndices.at(dih.idx1));
        auto &e_j = data->entry_at(particleIndices.at(dih.idx2));
        auto &e_k = data->entry_at(particleIndices.at(dih.idx3));
        auto &e_l = data->entry_at(particleIndices.at(dih.idx4));
        const auto x_ji = bcs::shortestDifference(e_j.pos, e_i.pos, context->boxSize(),
                                                  context->periodicBoundaryConditions());
        const auto x_kj = bcs::shortestDifference(e_k.pos, e_j.pos, context->boxSize(),
                                                  context->periodicBoundaryConditions());
        const auto x_kl = bcs::shortestDifference(e_k.pos, e_l.pos, context->boxSize(),
                                                  context->periodicBoundaryConditions());
        energy += potential->calculateEnergy(x_ji, x_kj, x_kl, dih);
        potential->calculateForce(e_i.force, e_j.force, e_k.force, e_l.force, x_ji, x_kj, x_kl, dih);
    }
    return energy;
}

namespace reactions {
namespace op {
CPUChangeParticleType::CPUChangeParticleType(CPUStateModel::data_type *const data,
                                             model::top::GraphTopology *const topology,
                                             const model::top::reactions::actions::TopologyReactionAction::vertex &v,
                                             const ParticleTypeId &type_to)
        : ChangeParticleType(topology, v, type_to), data(data) {}

void CPUChangeParticleType::execute() {
    const auto idx = topology->getParticles().at(_vertex->particleIndex);
    _vertex->setParticleType(previous_type);
    std::swap(data->entry_at(idx).type, previous_type);
}

void CPUChangeParticleType::undo() {
    execute();
}

}
}

}
}
}
}
}

