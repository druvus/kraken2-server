#pragma once

// A fixed-size pool of worker threads executing queued tasks. Built on
// ThreadSafeQueue; the destructor lets queued tasks finish, then joins.

#include <functional>
#include <thread>
#include <vector>

#include "thread_safe_queue.h"

class WorkerPool
{
public:
    // n == 0 selects std::thread::hardware_concurrency() (at least 1).
    explicit WorkerPool(size_t n = 0)
    {
        if (n == 0)
        {
            n = std::thread::hardware_concurrency();
            if (n == 0) n = 1;
        }
        for (size_t i = 0; i < n; i++)
        {
            threads_.emplace_back([this] {
                while (std::optional<std::function<void()>> task = tasks_.pop_wait())
                {
                    (*task)();
                }
            });
        }
    }

    ~WorkerPool()
    {
        tasks_.close();
        for (auto &t : threads_) t.join();
    }

    WorkerPool(const WorkerPool &) = delete;
    WorkerPool &operator=(const WorkerPool &) = delete;

    size_t thread_count() const { return threads_.size(); }

    // Queue a task. Tasks must not throw; wrap anything that might.
    void push(std::function<void()> task)
    {
        tasks_.push(std::move(task));
    }

    size_t queued() const { return tasks_.size(); }

private:
    ThreadSafeQueue<std::function<void()>> tasks_;
    std::vector<std::thread> threads_;
};
