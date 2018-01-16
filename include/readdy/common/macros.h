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


//
// Created by mho on 10/10/2016.
//

#pragma once

/**
 * Defines for the current OS
 */
#ifdef _WIN32
#define READDY_WINDOWS true
#define READDY_OSX false
#define READDY_LINUX false
#elif __APPLE__
#define READDY_WINDOWS false
#define READDY_OSX true
#define READDY_LINUX false
#elif __linux__ || __unix__
#define READDY_WINDOWS false
#define READDY_OSX false
#define READDY_LINUX true
#else
#error "Unknown compiler"
#endif

/**
 * Export symbols / change their visibility
 */
#if READDY_WINDOWS
#  define READDY_API __declspec(dllexport)
#else
#  define READDY_API __attribute__ ((visibility("default")))
#endif

#if !defined(NAMESPACE_BEGIN)
#  define NAMESPACE_BEGIN(name) namespace name {
#endif
#if !defined(NAMESPACE_END)
#  define NAMESPACE_END(name) }
#endif
