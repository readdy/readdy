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
 * Base class for all plugins.
 *
 * @file Plugin.h
 * @brief Header file containing the definitions for readdy::model::Plugin.
 * @author clonker
 * @date 04.03.16
 */

#pragma once
#include <string>
#include <memory>
#include <type_traits>
#include <sstream>
#include <unordered_map>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

class Plugin {
public:
    virtual const std::string &getName() const = 0;

    virtual ~Plugin() {};
};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
