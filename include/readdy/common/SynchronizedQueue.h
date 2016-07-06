/**
 * This queue synchronizes access.
 *
 * @file SynchronizedQueue.h
 * @brief Queue wrapper that synchronizes the access.
 * @author clonker
 * @date 28.06.16
 */

#ifndef READDY_MAIN_SYNCHRONIZEDQUEUE_H
#define READDY_MAIN_SYNCHRONIZEDQUEUE_H

#include <queue>
#include <mutex>
#include <condition_variable>

namespace readdy {
    namespace utils {
        template<typename E>
        class SynchronizedQueue {
        public:
            SynchronizedQueue() : queue(), mutex(), conditionVariable() {}

            void add(E e) {
                std::lock_guard<std::mutex> lock(mutex);
                queue.push(e);
                conditionVariable.notify_one();
            }

            E poll() {
                std::unique_lock<std::mutex> lock(mutex);
                while(queue.empty()) {
                    conditionVariable.wait(lock);
                }
                E e = queue.front();
                queue.pop();
                return e;
            }

        protected:
            std::queue<E> queue;
            mutable std::mutex mutex;
            std::condition_variable conditionVariable;
        };
    }
}

#endif //READDY_MAIN_SYNCHRONIZEDQUEUE_H
