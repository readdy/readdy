/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/

/**
 * @file Timer.h
 * @brief RAII timer with static storage and ability to dump results into a json string.
 * @author chrisfroe
 * @author clonker
 * @date 26.07.19
 */

#pragma once

#include <mutex>
#include <unordered_map>

namespace readdy::util {

struct PerformanceData {
    using time = double;
    /**
     * Creates a new performance datum containing the cumulative time and the number of calls.
     * @param t the initial cumulative time
     * @param c the initial number of calls
     */
    explicit PerformanceData(time t = 0, std::size_t c = 0) : _cumulativeTime(t), _count(c) {}

    /**
     * Record some elapsed time
     * @param elapsed the elapsed time
     */
    void record(time elapsed) const {
        std::unique_lock<std::mutex> lock(mutex);
        _cumulativeTime += elapsed;
        ++_count;
    }

    /**
     * the current value of the aggregated time in this datum
     * @return the cumulative time
     */
    time cumulativeTime() const {
        return _cumulativeTime;
    }

    /**
     * the number of calls that were made to record plus the initial number of counts given in the constructor
     * @return the number of calls
     */
    std::size_t count() const {
        return _count;
    }

    /**
     * clears this datum
     */
    void clear() const {
        std::unique_lock<std::mutex> lock(mutex);
        _cumulativeTime = 0.;
        _count = 0;
    }

private:
    mutable time _cumulativeTime = 0.;
    mutable std::size_t _count = 0;
    mutable std::mutex mutex;
};

class Timer {
public:
    /**
     * Creates a new timer object that will record the time elapsed between construction and destruction in a given
     * performance datum.
     *
     * @param target the performance datum to store the elapsed time in
     * @param measure if set to false, the measurements will not be stored
     */
    Timer(const PerformanceData &target, bool measure) : target(target), measure(measure) {
        if (measure) {
            begin = std::chrono::high_resolution_clock::now();
        }
    }

    explicit Timer(const std::string &targetName, bool measure = true) : Timer(perf[targetName], measure) {};

    ~Timer() {
        stop();
    }

    Timer(const Timer &) = delete;

    Timer(Timer &&) = delete;

    Timer &operator=(const Timer &) = delete;

    Timer &operator=(Timer &&) = delete;

    /**
     * Recording the elapsed time since construction of this object into the given performance datum.
     * In cases where scoped timing is not appropriate stop() can be called manually, otherwise it is called
     * in the destructor. */
    void stop() {
        if (measure && !wasMeasured) {
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            long elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
            auto elapsedSeconds =
                    static_cast<PerformanceData::time>(1e-6) * static_cast<PerformanceData::time>(elapsed);
            target.record(elapsedSeconds);
            wasMeasured = true;
        }
    }

    static std::string perfToJsonString();

    static void clear();

private:
    bool measure{true};
    bool wasMeasured{false};
    const PerformanceData &target;
    std::chrono::high_resolution_clock::time_point begin;
    static std::unordered_map<std::string, PerformanceData> perf;
};

}
