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
 * @file TrajectoryObservable.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.08.16
 */

#ifndef READDY_MAIN_TRAJECTORYOBSERVABLE_H
#define READDY_MAIN_TRAJECTORYOBSERVABLE_H

#include <readdy/model/observables/Observables.h>

namespace readdy {
    namespace io {
        class TrajectoryWriter : public readdy::model::observables::Combiner<void, readdy::model::observables::ParticlePosition> {
            using kernel_t = readdy::model::Kernel;
            using ppObs_t = readdy::model::observables::ParticlePosition;
        public:
            TrajectoryWriter(const std::string& path, kernel_t *const kernel, unsigned int stride, ppObs_t* parent)
                    : readdy::model::observables::Combiner(kernel, stride, parent) {};

            virtual void evaluate() = 0;

        protected:

        };
    }
}
#endif //READDY_MAIN_TRAJECTORYOBSERVABLE_H
