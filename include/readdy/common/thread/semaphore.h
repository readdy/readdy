/**
 * << detailed description >>
 *
 * @file semaphore.h
 * @brief << brief description >>
 * @author clonker
 * @date 16.11.16
 */

#ifndef READDY_MAIN_SEMAPHORE_H
#define READDY_MAIN_SEMAPHORE_H

#include <mutex>
#include <condition_variable>


namespace readdy {
namespace util {
namespace thread {

/**
 * counting semaphore implementation
 */
class semaphore {
public:
    explicit semaphore(int count = 0) : count(count) {}

    /**
     * increase count
     */
    inline void signal() const {
        std::unique_lock<std::mutex> lock(mutex);
        count++;
        cv.notify_one();
    }

    /**
     * waste time while count is <= 0, then decrement
     */
    inline void wait() const {
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock, [this]() { return count > 0; });
        count--;
    }

    /**
     * waste time while count is <= 0
     */
    inline void wait_no_decrement() const {
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock, [this]() { return count > 0; });
    }

    semaphore(const semaphore&) = delete;
    semaphore(semaphore&&) = delete;
    semaphore& operator=(const semaphore&) = delete;
    semaphore& operator=(semaphore&&) = delete;

private:
    mutable std::mutex mutex;
    mutable std::condition_variable cv;
    mutable int count;
};

}
}
}
#endif //READDY_MAIN_SEMAPHORE_H
