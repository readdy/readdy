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
 * Including this header will expose the C++ top-level API, i.e., everything under
 *   * ObservableHandle.h
 *   * Simulation.h
 *   * SimulationScheme.h
 *
 * @file readdy.h
 * @brief Main header of ReaDDy. Including this will expose the C++ top-level API.
 * @author clonker
 * @date 06.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "api/ObservableHandle.h"
#include "api/Simulation.h"
#include "readdy/api/SimulationLoop.h"
