/*********************************************************
*
*  Copyright (C) 2014 by Vitaliy Vitsentiy
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*********************************************************/


#pragma once

#include <functional>
#include <thread>
#include <atomic>
#include <vector>
#include <memory>
#include <exception>
#include <future>
#include <mutex>
#include <queue>

// thread pool to run user's functors with signature
//      ret func(int id, other_params)
// where id is the index of the thread that runs the functor
// ret is some return type

namespace ctpl {

namespace detail {
template<typename T>
class Queue {
public:

    template<typename Iterator>
    bool pushAll(Iterator begin, Iterator end) {
        std::unique_lock<std::mutex> lock(this->mutex);
        std::for_each(begin, end, [this](const auto &val) {
            this->q.push(val);
        });
        return true;
    }

    bool push(T const &value) {
        std::unique_lock<std::mutex> lock(this->mutex);
        this->q.push(value);
        return true;
    }

    // deletes the retrieved element, do not use for non integral types
    bool pop(T &v) {
        std::unique_lock<std::mutex> lock(this->mutex);
        if (this->q.empty())
            return false;
        v = this->q.front();
        this->q.pop();
        return true;
    }

    void popAll() {
        std::unique_lock<std::mutex> lock(this->mutex);
        while(!q.empty()) q.pop();
    }

    bool empty() {
        std::unique_lock<std::mutex> lock(this->mutex);
        return this->q.empty();
    }

    Queue() = default;

private:
    std::queue<T> q;
    std::mutex mutex;
};
}

class thread_pool {

public:

    thread_pool() : thread_pool(0) {  }

    explicit thread_pool(int nThreads) : nWaiting(0), isStop(false), isDone(false) {
        this->resize(static_cast<std::size_t>(nThreads));
    }

    // the destructor waits for all the functions in the queue to be finished
    ~thread_pool() {
        this->stop(true);
    }

    // get the number of running threads in the pool
    std::size_t size() const { return this->threads.size(); }

    // number of idle threads
    int n_idle() { return this->nWaiting; }

    std::thread &get_thread(int i) { return *this->threads[i]; }

    void resize_wait(std::size_t n) {
        stop(true);
        nWaiting = 0;
        isStop = false;
        isDone = false;
        resize(n);
    }

    // change the number of threads in the pool
    // should be called from one thread, otherwise be careful to not interleave, also with this->stop()
    // nThreads must be >= 0
    void resize(std::size_t nThreads) {
        if (!this->isStop && !this->isDone) {
            auto oldNThreads = static_cast<int>(this->threads.size());
            if (oldNThreads <= nThreads) {  // if the number of threads is increased
                this->threads.resize(nThreads);
                this->flags.resize(nThreads);

                for (int i = oldNThreads; i < nThreads; ++i) {
                    this->flags[i] = std::make_shared<std::atomic<bool>>(false);
                    this->set_thread(i);
                }
            } else {  // the number of threads is decreased
                for (int i = oldNThreads - 1; i >= nThreads; --i) {
                    *this->flags[i] = true;  // this thread will finish
                    this->threads[i]->detach();
                }
                {
                    // stop the detached threads that were waiting
                    std::unique_lock<std::mutex> lock(this->mutex);
                    this->cv.notify_all();
                }
                // safe to delete because the threads are detached
                this->threads.resize(nThreads);
                // safe to delete because the threads have copies of shared_ptr of the flags, not originals
                this->flags.resize(nThreads);
            }
        }
    }

    // empty the queue
    void clear_queue() {
        std::function<void(int id)> *_f;
        while (this->q.pop(_f))
            delete _f; // empty the queue
    }

    // pops a functional wrapper to the original function
    std::function<void(int)> pop() {
        std::function<void(int id)> *_f = nullptr;
        this->q.pop(_f);
        // at return, delete the function even if an exception occurred
        std::unique_ptr<std::function<void(int id)>> func(_f);
        std::function<void(int)> f;
        if (_f) f = *_f;
        return f;
    }

    template<typename F, typename... Args>
    auto pack(F &&f, Args &&... args) const -> std::function<decltype(f(0, args...))(std::size_t)> {
        return std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Args>(args)...);
    }

    // wait for all computing threads to finish and stop all threads
    // may be called asynchronously to not pause the calling thread while waiting
    // if isWait == true, all the functions in the queue are run, otherwise the queue is cleared without running the functions
    void stop(bool isWait = false) {
        if (!isWait) {
            if (this->isStop)
                return;
            this->isStop = true;
            for (int i = 0, n = static_cast<int>(this->size()); i < n; ++i) {
                *this->flags[i] = true;  // command the threads to stop
            }
            this->clear_queue();  // empty the queue
        } else {
            if (this->isDone || this->isStop)
                return;
            this->isDone = true;  // give the waiting threads a command to finish
        }
        {
            std::unique_lock<std::mutex> lock(this->mutex);
            this->cv.notify_all();  // stop all waiting threads
        }
        for (auto &thread : this->threads) {
            // wait for the computing threads to finish
            if (thread->joinable()) {
                thread->join();
            }
        }
        // if there were no threads in the pool but some functors in the queue, the functors are not deleted by the threads
        // therefore delete them here
        this->clear_queue();
        this->threads.clear();
        this->flags.clear();
    }

    template<typename F>
    auto pushAll(std::vector<F> &&funs) -> std::vector<std::future<decltype(std::declval<F>()(0))>> {
        using ResultVec = std::vector<std::future<decltype(std::declval<F>()(0))>>;
        ResultVec vec;
        vec.reserve(funs.size());

        std::vector<std::function<void(int)>*> internalFunctions;
        internalFunctions.reserve(funs.size());

        for(auto &&f : funs) {
            auto pck = std::make_shared<std::packaged_task<decltype(f(0))(int)>>(std::forward<F>(f));
            auto _f = new std::function<void(int id)>([pck](int id) {
                (*pck)(id);
            });
            internalFunctions.push_back(_f);
            vec.push_back(pck->get_future());
        }

        this->q.pushAll(internalFunctions.begin(), internalFunctions.end());
        std::unique_lock<std::mutex> lock(this->mutex);
        for(auto i=0U; i < funs.size(); ++i) this->cv.notify_one();
        return vec;
    }

    template<typename F, typename... Rest>
    auto push(F &&f, Rest &&... rest) -> std::future<decltype(f(0, rest...))> {
        auto pck = std::make_shared<std::packaged_task<decltype(f(0, rest...))(int)>>(
                std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Rest>(rest)...)
        );
        auto _f = new std::function<void(int id)>([pck](int id) {
            (*pck)(id);
        });
        this->q.push(_f);
        std::unique_lock<std::mutex> lock(this->mutex);
        this->cv.notify_one();
        return pck->get_future();
    }

    // run the user's function that excepts argument int - id of the running thread. returned value is templatized
    // operator returns std::future, where the user can get the result and rethrow the catched exceptins
    template<typename F>
    auto push(F &&f) -> std::future<decltype(f(0))> {
        auto pck = std::make_shared<std::packaged_task<decltype(f(0))(int)>>(std::forward<F>(f));
        auto _f = new std::function<void(int id)>([pck](int id) {
            (*pck)(id);
        });
        this->q.push(_f);
        std::unique_lock<std::mutex> lock(this->mutex);
        this->cv.notify_one();
        return pck->get_future();
    }

    thread_pool(const thread_pool &) = delete;// = delete;
    thread_pool(thread_pool &&) = delete;
    thread_pool &operator=(const thread_pool &) = delete;// = delete;
    thread_pool &operator=(thread_pool &&) = delete;
private:


    void set_thread(int i) {
        std::shared_ptr<std::atomic<bool>> flag(this->flags[i]); // a copy of the shared ptr to the flag
        auto f = [this, i, flag/* a copy of the shared ptr to the flag */]() -> void {
            std::atomic<bool> &_flag = *flag;
            std::function<void(int id)> *_f;
            bool isPop = this->q.pop(_f);
            while (true) {
                while (isPop) {  // if there is anything in the queue
                    std::unique_ptr<std::function<void(int id)>> func(
                            _f); // at return, delete the function even if an exception occurred
                    (*_f)(i);
                    if (_flag)
                        return;  // the thread is wanted to stop, return even if the queue is not empty yet
                    else
                        isPop = this->q.pop(_f);
                }
                // the queue is empty here, wait for the next command
                std::unique_lock<std::mutex> lock(this->mutex);
                ++this->nWaiting;
                this->cv.wait(lock, [this, &_f, &isPop, &_flag]() {
                    isPop = this->q.pop(_f);
                    return isPop || this->isDone || _flag;
                });
                --this->nWaiting;
                if (!isPop)
                    return;  // if the queue is empty and this->isDone == true or *flag then return
            }
        };
        this->threads[i] = std::make_unique<std::thread>(f);
    }

    std::vector<std::unique_ptr<std::thread>> threads;
    std::vector<std::shared_ptr<std::atomic<bool>>> flags;
    detail::Queue<std::function<void(int id)> *> q;
    std::atomic<bool> isDone;
    std::atomic<bool> isStop;
    std::atomic<int> nWaiting;  // how many threads are waiting

    std::mutex mutex;
    std::condition_variable cv;
};

}
