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

namespace readdy {
namespace util {
struct Timer {

    Timer(std::string label, bool print = true) : label(label), print(print) {}

    ~Timer() {
        if (print) {
            std::cout << "Elapsed (" << label << "): " << getSeconds() << " seconds" << std::endl;
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
