/**
 * << detailed description >>
 *
 * @file ScopedThread.h
 * @brief << brief description >>
 * @author clonker
 * @date 01.08.16
 */

#ifndef READDY_CPUKERNEL_SCOPEDTHREAD_H
#define READDY_CPUKERNEL_SCOPEDTHREAD_H

#include <thread>

namespace readdy {
namespace kernel {
namespace cpu {
namespace util {
class scoped_thread {
    std::thread t;
public:
    explicit scoped_thread(std::thread _t) : t(std::move(_t)) {
        if (!t.joinable()) throw std::logic_error("No thread!");
    }

    ~scoped_thread() {
        if (t.joinable()) {
            t.join();
        }
    }

    scoped_thread(scoped_thread &&rhs) {
        t = std::move(rhs.t);
    }

    scoped_thread &operator=(scoped_thread &&rhs) {
        t = std::move(rhs.t);
        return *this;
    }

    scoped_thread(const scoped_thread &) = delete;

    scoped_thread &operator=(const scoped_thread &) = delete;
};
}
}
}
}
#endif //READDY_CPUKERNEL_SCOPEDTHREAD_H
