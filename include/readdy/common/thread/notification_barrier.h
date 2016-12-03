/**
 * << detailed description >>
 *
 * @file notification_barrier.h
 * @brief << brief description >>
 * @author clonker
 * @date 24.11.16
 */

#ifndef READDY_MAIN_NOTIFICATION_BARRIER_H
#define READDY_MAIN_NOTIFICATION_BARRIER_H

#include <mutex>
#include <condition_variable>
#include <atomic>

namespace readdy {
namespace util {
namespace thread {

class notification_barrier {
public:
    explicit notification_barrier() { }

    /**
     * waste time until all threads reached the barrier
     */
    void wait() const {
        std::unique_lock<std::mutex> lock(mutex);
        while(!done.load()) {
            cv.wait(lock, [this] { return done.load(); });
        }
    }

    void ready() const {
        done.store(true);
        cv.notify_all();
    }

private:
    mutable std::mutex mutex;
    mutable std::condition_variable cv;
    mutable std::atomic<bool> done {false};
};

}
}
}
#endif //READDY_MAIN_NOTIFICATION_BARRIER_H
