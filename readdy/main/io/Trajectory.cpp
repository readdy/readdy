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
 * @file Trajectory.cpp
 * @brief Core library impl for Trajectory.h
 * @author chrisfroe
 * @date 12.12.16
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/io/Trajectory.h>

namespace obs = readdy::model::observables;

namespace readdy {
namespace io {

Trajectory::Trajectory(model::Kernel *const kernel, unsigned int stride, unsigned int flushStride, obs::Particles *particlesObservable)
        : obs::Combiner<std::pair<std::vector<model::observables::time_step_type>, std::vector<obs::Particles::result_t>>, obs::Particles>(kernel,
                                                                                                                                           stride,
                                                                                                                                           particlesObservable) {};

void Trajectory::evaluate() {
    const auto &currentInput = std::get<0>(this->parentObservables)->getResult();
    auto &resultTimes = std::get<0>(this->result);
    auto &resultValues = std::get<1>(this->result);
    resultTimes.push_back(this->getCurrentTimeStep());
    resultValues.push_back(currentInput);
    if (count % flushStride == 0) {
        //flush();
    }
    ++count;
}

void Trajectory::flush() {
//writer.write(this->result);
    auto &resultTimes = std::get<0>(this->result);
    auto &resultValues = std::get<1>(this->result);
    resultTimes.clear();
    resultValues.clear();
}

}
}