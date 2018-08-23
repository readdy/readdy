/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
        if(bond.forceConstant == 0) continue;

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
    _vertex->particleType() = previous_type;
    std::swap(data->entry_at(idx).type, previous_type);
}

void CPUChangeParticleType::undo() {
    execute();
}

CPUChangeParticlePosition::CPUChangeParticlePosition(
        CPUStateModel::data_type *const data, model::top::GraphTopology *topology,
        const model::top::reactions::actions::TopologyReactionAction::vertex &v, Vec3 position)
        : ChangeParticlePosition(topology, v, position), data(data) { }

void CPUChangeParticlePosition::execute() {
    const auto idx = topology->getParticles().at(_vertex->particleIndex);
    std::swap(data->entry_at(idx).pos, _posTo);
}

void CPUChangeParticlePosition::undo() {
    execute();
}

}
}

}
}
}
}
}

