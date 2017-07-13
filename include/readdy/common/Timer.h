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
 * Header file for the Timer class, which measures execution time in RAII style.
 *
 * @file Timer.h
 * @brief Header file containing definition and implementation of a RAII-Timer.
 * @author clonker
 * @date 13.07.16
 */

#pragma once
#include <memory>
#include <chrono>

#include "logging.h"
#include "common.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
/**
 * Timer class.
 */
class Timer {
public:
    /**
     * constructs a new timer
     * @param label the label of the timer
     * @param print if the timer should print the elapsed time when it runs out of scope
     */
    Timer(const std::string &label, bool print = true) : label(label), print(print) {}

    /**
     * the destructor of the Timer class. prints the elapsed time if Timer::print is true.
     */
    ~Timer() {
        if (print) {
            log::debug("Elapsed ({}): {} seconds", label, getSeconds());
        }
    }

    /**
     * returns the elapsed seconds
     * @return the elapsed seconds
     */
    scalar getSeconds() {
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        long elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        return (scalar) 1e-6 * (scalar) elapsed;
    }

private :
    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
    std::string label;
    bool print;
};
NAMESPACE_END(util)
NAMESPACE_END(readdy)
