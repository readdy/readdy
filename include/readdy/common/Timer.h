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
 * @author chrisfroe
 * @date 13.07.16
 */

#pragma once

#include <memory>
#include <chrono>

#include "logging.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

class Timer {
public:
    using cumulative_time_map = std::unordered_map<std::string, double>;
    using counts_map = std::unordered_map<std::string, std::size_t>;

    Timer(const std::string &label) : label(label) {}

    void start() {
        begin = std::chrono::high_resolution_clock::now();
    }

    void stop() const {
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        long elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        const auto elapsedSeconds = (double) 1e-6 * (double) elapsed;
        std::unique_lock<std::mutex> lock(mutex);
        cumulativeTime[label] += elapsedSeconds;
        _counts[label]++;
    }

    static const cumulative_time_map &times() {
        return cumulativeTime;
    }

    static const counts_map &counts() {
        return _counts;
    }

private:
    std::string label;
    std::chrono::high_resolution_clock::time_point begin;
    static cumulative_time_map cumulativeTime;
    static counts_map _counts;
    static std::mutex mutex;
};

class RAIITimer {
public:
    /**
     * constructs a new timer
     * @param label the label of the timer
     * @param measure if the timer should measure, if false, it does nothing
     */
    RAIITimer(bool measure, const std::string &label) : _measure(measure), timer(Timer(label)) {
        if (_measure) {
            timer.start();
        }
    }

    ~RAIITimer() {
        if (_measure) {
            timer.stop();
        }
    }

private:
    Timer timer;
    bool _measure;
};
NAMESPACE_END(util)
NAMESPACE_END(readdy)
