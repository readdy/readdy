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

#ifndef READDY_MAIN_TIMER_H
#define READDY_MAIN_TIMER_H

#include <memory>
#include <chrono>

#include "logging.h"

namespace readdy {
namespace util {
struct Timer {

    Timer(const std::string &label, bool print = true) : label(label), print(print) {}

    ~Timer() {
        if (print) {
            log::debug("Elapsed ({}): {} seconds", label, getSeconds());
        }
    }

    double getSeconds() {
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        long elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        return (double) 1e-6 * (double) elapsed;
    }

private :
    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
    std::string label;
    bool print;
};
}
}
#endif //READDY_MAIN_TIMER_H
