/**
 * << detailed description >>
 *
 * @file ScopedThread.h
 * @brief << brief description >>
 * @author clonker
 * @date 01.08.16
 */

#ifndef READDY_MAIN_SCOPEDTHREAD_H
#define READDY_MAIN_SCOPEDTHREAD_H

#include <thread>
#include <boost/log/trivial.hpp>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace util {
                class ScopedThread {
                    std::thread t;
                public:
                    explicit ScopedThread(std::thread _t) : t(std::move(_t)) {
                        if (!t.joinable()) throw std::logic_error("No thread!");
                    }

                    ~ScopedThread() {
                        if(t.joinable()) t.join();
                    }

                    ScopedThread(ScopedThread&& rhs) {
                        t = std::move(rhs.t);
                    }
                    ScopedThread& operator=(ScopedThread&& rhs) {
                        t = std::move(rhs.t);
                        return *this;
                    }

                    ScopedThread(const ScopedThread &) = delete;
                    ScopedThread &operator=(const ScopedThread &) = delete;
                };
            }
        }
    }
}
#endif //READDY_MAIN_SCOPEDTHREAD_H
