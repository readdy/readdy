/**
 * << detailed description >>
 *
 * @file ThreadGuard.h
 * @brief << brief description >>
 * @author clonker
 * @date 01.08.16
 */

#ifndef READDY_CPUKERNEL_THREADGUARD_H
#define READDY_CPUKERNEL_THREADGUARD_H

#include <thread>

namespace readdy {
namespace kernel {
namespace cpu {
namespace util {
class ThreadGuard {
    std::thread *t;
public:
    explicit ThreadGuard(std::thread *const t) : t(t) {}

    ~ThreadGuard() {
        if (t && t->joinable()) t->join();
    }

    ThreadGuard(const ThreadGuard &) = delete;

    ThreadGuard &operator=(const ThreadGuard &) = delete;

    ThreadGuard(ThreadGuard &&tg) : t(std::move(tg.t)) {
        tg.t = nullptr;
    };

    ThreadGuard &operator=(ThreadGuard &&tg) {
        t = std::move(tg.t);
        tg.t = nullptr;
        return *this;
    }
};
}
}
}
}
#endif //READDY_CPUKERNEL_THREADGUARD_H
