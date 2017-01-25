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
 * @file Utils.h
 * @brief << brief description >>
 * @author clonker
 * @date 17.11.16
 */

#ifndef READDY_MAIN_MODEL_UTILS_H
#define READDY_MAIN_MODEL_UTILS_H

#include <readdy/model/KernelContext.h>
#include <readdy/model/observables/Observable.h>

namespace readdy {
namespace model {
namespace util {

double getRecommendedTimeStep(unsigned int N, KernelContext&);
double getMaximumDisplacement(KernelContext&, const double);

}
}
}

#endif //READDY_MAIN_MODEL_UTILS_H
