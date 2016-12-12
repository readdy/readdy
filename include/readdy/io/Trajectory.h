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
 * @file Trajectory.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 12.12.16
 * @copyright GNU Lesser General Public License v3.0
 */
#ifndef READDY_MAIN_TRAJECTORY_H
#define READDY_MAIN_TRAJECTORY_H

#include <array>
#include <readdy/model/Kernel.h>

namespace readdy {
namespace io {

class Trajectory
        : public model::observables::Combiner<std::pair<std::vector<model::observables::time_step_type>, std::vector<model::observables::Particles::result_t>>, model::observables::Particles> {

public:
    Trajectory(model::Kernel *const kernel, unsigned int stride, unsigned int flushStride, model::observables::Particles *particlesObservable);

    virtual void evaluate();

    virtual void flush();

protected:
    unsigned int flushStride;
    unsigned int count = 0;
};

}
}

#endif //READDY_MAIN_TRAJECTORY_H
