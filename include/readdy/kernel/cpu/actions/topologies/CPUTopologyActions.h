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
 * @file CPUTopologyActions.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.02.17
 * @copyright BSD-3
 */

#pragma once

#include <readdy/model/topologies/Topology.h>
#include <readdy/model/topologies/potentials/TopologyPotentialActions.h>
#include <readdy/kernel/cpu/CPUStateModel.h>
#include <readdy/model/topologies/GraphTopology.h>

namespace readdy::kernel::cpu::actions::top {

class CPUCalculateHarmonicBondPotential : public readdy::model::top::pot::CalculateHarmonicBondPotential {

    const harmonic_bond *const potential;
    CPUStateModel::data_type *const data;

public:
    CPUCalculateHarmonicBondPotential(const readdy::model::Context *context,
                                      CPUStateModel::data_type *data,
                                      const harmonic_bond *potential);

    readdy::scalar perform(const readdy::model::top::GraphTopology *topology) override;
};


class CPUCalculateHarmonicAnglePotential : public readdy::model::top::pot::CalculateHarmonicAnglePotential {
    const harmonic_angle *const potential;
    CPUStateModel::data_type *const data;
public:
    CPUCalculateHarmonicAnglePotential(const readdy::model::Context *context, CPUStateModel::data_type *data,
                                       const harmonic_angle *potential);

    readdy::scalar perform(const readdy::model::top::GraphTopology* topology) override;
};

class CPUCalculateCosineDihedralPotential : public readdy::model::top::pot::CalculateCosineDihedralPotential {
    const cos_dihedral *const potential;
    CPUStateModel::data_type *const data;
public:
    CPUCalculateCosineDihedralPotential(const readdy::model::Context *context,
                                        CPUStateModel::data_type *data,
                                        const cos_dihedral *pot);

    readdy::scalar perform(const readdy::model::top::GraphTopology* topology) override;
};

namespace reactions::op {

class CPUChangeParticleType : public readdy::model::top::reactions::actions::ChangeParticleType {
    CPUStateModel::data_type *const data;
public:
    CPUChangeParticleType(CPUStateModel::data_type *data, model::top::GraphTopology *topology,
                          const readdy::model::top::Graph::PersistentVertexIndex &v,
                          const ParticleTypeId &type_to);

    void execute() override;

};

class CPUChangeParticlePosition : public readdy::model::top::reactions::actions::ChangeParticlePosition {
    CPUStateModel::data_type *const data;
public:
    CPUChangeParticlePosition(CPUStateModel::data_type *data, model::top::GraphTopology *topology,
                              const readdy::model::top::Graph::PersistentVertexIndex &v, Vec3 position);

    void execute() override;
};

class CPUAppendParticle : public readdy::model::top::reactions::actions::AppendParticle {
    CPUStateModel::data_type *const data;
    readdy::model::Particle particle;
    std::size_t insertIndex {0};
public:
    CPUAppendParticle(CPUStateModel::data_type *const data, model::top::GraphTopology *topology,
                      std::vector<readdy::model::top::Graph::PersistentVertexIndex> neighbors, ParticleTypeId type, Vec3 pos)
            : AppendParticle(topology, std::move(neighbors), type, pos), data(data), particle(pos, type) {};

    void execute() override {
        auto entry = CPUStateModel::data_type::entry_type(particle);
        insertIndex = data->addEntry(std::move(entry));

        auto firstNeighbor = neighbors[0];
        auto ix = topology->appendParticle(insertIndex, firstNeighbor);
        for(auto v : neighbors) {
            topology->addEdge(v, ix);
        }
    }
};

}

}
