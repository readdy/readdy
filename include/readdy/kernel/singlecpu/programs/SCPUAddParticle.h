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
 * @file SingleCPUAddParticleProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_SINGLECPUADDPARTICLEPROGRAM_H
#define READDY_MAIN_SINGLECPUADDPARTICLEPROGRAM_H

#include <readdy/model/programs/Programs.h>
#include <readdy/model/Particle.h>
#include <vector>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {
class SCPUAddParticle : public readdy::model::programs::AddParticle {
public:
    SCPUAddParticle(SCPUKernel *kernel);

    virtual ~SCPUAddParticle() override;

    virtual void execute() override;

    void addParticle(readdy::model::Particle particle) {
        particles.emplace_back(std::move(particle));
    }

    void setParticles(std::vector<readdy::model::Particle> particles) {
        this->particles = std::move(particles);
    }

private:
    std::vector<readdy::model::Particle> particles;
    SCPUKernel *kernel;
};
}
}
}
}


#endif //READDY_MAIN_SINGLECPUADDPARTICLEPROGRAM_H
