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
 * Nodelete callable struct that can be supplied to, e.g., instances of unique_ptr such that the helt object is not
 * deleted along with the instance.
 *
 * @file nodelete.h
 * @brief Header file containing nodelete.
 * @author clonker
 * @date 18.10.16
 */

#pragma once
NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
struct nodelete {
    template<typename T>
    void operator()(T *) {}
};
NAMESPACE_END(util)
NAMESPACE_END(readdy)
