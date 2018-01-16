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
 * Header file containing definitions for various observables. Currently:
 *  - Positions,
 *  - Particles,
 *  - RadialDistribution,
 *  - HistogramAlongAxis,
 *  - NParticles,
 *  - Forces,
 *  - Reactions,
 *  - ReactionCounts
 *
 * @file Observables.h
 * @brief Header file combining definitions for various observables.
 * @author clonker
 * @date 26.04.16
 * @todo for demonstration purposes, add a more meaningful combiner observable, such as velocity
 */

#pragma once

#include "HistogramAlongAxis.h"
#include "Particles.h"
#include "Positions.h"
#include "RadialDistribution.h"
#include "Forces.h"
#include "NParticles.h"
#include "Reactions.h"
#include "ReactionCounts.h"
#include "Energy.h"
#include "io/Trajectory.h"
