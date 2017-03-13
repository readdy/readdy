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
 * @file SCPUTopologyActions.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include <readdy/model/topologies/actions/TopologyActions.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>
#include <readdy/model/topologies/Topology.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(kernel)
NAMESPACE_BEGIN(scpu)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class SCPUCalculateHarmonicBondPotential : public readdy::model::top::CalculateHarmonicBondPotential {

    const readdy::model::top::HarmonicBondPotential *const potential;
    SCPUParticleData *const data;

public:
    SCPUCalculateHarmonicBondPotential(const readdy::model::KernelContext *const context,
                                       SCPUParticleData *const data,
                                       const readdy::model::top::HarmonicBondPotential *const potential)
            : CalculateHarmonicBondPotential(context), potential(potential), data(data) {}

    virtual double perform() override {
        readdy::model::Vec3::value_t energy = 0;
        const auto& particleIndices = potential->getTopology()->getParticles();
        const auto& d = context->getShortestDifferenceFun();
        for(const auto& bond : potential->getBonds()) {
            readdy::model::Vec3 forceUpdate{0, 0, 0};
            auto& e1 = data->entry_at(particleIndices.at(bond.idx1));
            auto& e2 = data->entry_at(particleIndices.at(bond.idx2));
            const auto x_ij = d(e1.position(), e2.position());
            potential->calculateForce(forceUpdate, x_ij, bond);
            e1.force += forceUpdate;
            e2.force += -1 * forceUpdate;
            energy += potential->calculateEnergy(x_ij, bond);
        }
        return energy;
    }

};

class SCPUCalculateHarmonicAnglePotential : public readdy::model::top::CalculateHarmonicAnglePotential {
    const readdy::model::top::HarmonicAnglePotential *const potential;
    SCPUParticleData *const data;
public:
    SCPUCalculateHarmonicAnglePotential(const readdy::model::KernelContext *const context, SCPUParticleData *const data,
                                        const readdy::model::top::HarmonicAnglePotential*const potential)
            : CalculateHarmonicAnglePotential(context), potential(potential), data(data) {}

    virtual double perform() override {
        readdy::model::Vec3::value_t energy = 0;
        const auto& particleIndices = potential->getTopology()->getParticles();
        const auto& d = context->getShortestDifferenceFun();


        for(const auto& angle : potential->getAngles()) {
            auto& e1 = data->entry_at(particleIndices.at(angle.idx1));
            auto& e2 = data->entry_at(particleIndices.at(angle.idx2));
            auto& e3 = data->entry_at(particleIndices.at(angle.idx3));
            const auto x_ji = d(e2.pos, e1.pos);
            const auto x_jk = d(e2.pos, e3.pos);
            energy += potential->calculateEnergy(x_ji, x_jk, angle);
            potential->calculateForce(e1.force, e2.force, e3.force, x_ji, x_jk, angle);
        }
        return energy;
    }

};

class SCPUCalculateCosineDihedralPotential : public readdy::model::top::CalculateCosineDihedralPotential {
    const readdy::model::top::CosineDihedralPotential *const potential;
    SCPUParticleData *const data;
public:
    SCPUCalculateCosineDihedralPotential(const readdy::model::KernelContext *const context,
                                         SCPUParticleData *const data,
                                         const readdy::model::top::CosineDihedralPotential* const pot)
            : CalculateCosineDihedralPotential(context), potential(pot), data(data){
    }

    virtual double perform() override {
        readdy::model::Vec3::value_t energy = 0;
        const auto& particleIndices = potential->getTopology()->getParticles();
        const auto& d = context->getShortestDifferenceFun();

        for(const auto& dih : potential->getDihedrals()) {
            auto& e_i = data->entry_at(particleIndices.at(dih.idx1));
            auto& e_j = data->entry_at(particleIndices.at(dih.idx2));
            auto& e_k = data->entry_at(particleIndices.at(dih.idx3));
            auto& e_l = data->entry_at(particleIndices.at(dih.idx4));
            const auto x_ji = d(e_j.pos, e_i.pos);
            const auto x_kj = d(e_k.pos, e_j.pos);
            const auto x_kl = d(e_k.pos, e_l.pos);
            energy += potential->calculateEnergy(x_ji, x_kj, x_kl, dih);
            potential->calculateForce(e_i.force, e_j.force, e_k.force, e_l.force, x_ji, x_kj, x_kl, dih);
        }
        return energy;
    }
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(scpu)
NAMESPACE_END(kernel)
NAMESPACE_END(readdy)
