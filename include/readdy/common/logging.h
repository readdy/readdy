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
 * Header file containing a console() logger which gives the spd console logger.
 *
 * @file logging.h
 * @brief Header responsible for the logging system.
 * @author clonker
 * @date 14.10.16
 */

#pragma once

#include <spdlog/spdlog.h>
#include "macros.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(log)

std::shared_ptr<spdlog::logger> console();

template<typename... Args>
void trace(Args &&... args) {
    console()->trace(std::forward<Args>(args)...);
}

template<typename... Args>
void debug(Args &&... args) {
    console()->debug(std::forward<Args>(args)...);
}

template<typename... Args>
void critical(Args &&... args) {
    console()->critical(std::forward<Args>(args)...);
}

template<typename... Args>
void warn(Args &&... args) {
    console()->warn(std::forward<Args>(args)...);
}

template<typename... Args>
void error(Args &&... args) {
    console()->error(std::forward<Args>(args)...);
}

template<typename... Args>
void info(Args &&... args) {
    console()->info(std::forward<Args>(args)...);
}

class Level {
public:
    Level(spdlog::level::level_enum newLevel);

    ~Level();

private:
    spdlog::level::level_enum oldLevel;
};

NAMESPACE_END(log)
NAMESPACE_END(readdy)
