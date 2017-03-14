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
 * @file Util_bits.h
 * @brief << brief description >>
 * @author clonker
 * @date 06.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/io/Group.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(io)
NAMESPACE_BEGIN(util)
inline bool groupExists(const Group &cwd, const std::string &name) {
    hid_t hid = cwd.getHandle();
    H5E_BEGIN_TRY
        {
            hid = H5Gopen(hid, name.c_str(), H5P_DEFAULT);
            if (hid > 0) {
                H5Gclose(hid);
            }
        }
    H5E_END_TRY
    return (hid > 0);
}

NAMESPACE_END(util)
NAMESPACE_END(io)
NAMESPACE_END(readdy)
