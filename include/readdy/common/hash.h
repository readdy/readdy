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
 * Defines utility functions for dealing with hashes.
 *
 * @file hash.h
 * @brief Utility functions for hashes.
 * @author clonker
 * @date 14.10.16
 */

#pragma once

#include <cstddef>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(hash)

/**
 * Simplified version of boost hash combine.
 * @param seed the seed
 * @param v the value
 */
template<typename T>
void combine(std::size_t &seed, const T &v) {
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct EnumClassHash {
    template<typename T>
    std::size_t operator()(T t) const {
        return static_cast<std::size_t>(t);
    }
};

NAMESPACE_END(hash)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
