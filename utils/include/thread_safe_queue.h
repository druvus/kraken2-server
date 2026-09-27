#pragma once

#include <condition_variable>
#include <limits>
#include <mutex>
#include <optional>
#include <queue>
#include <utility>

// A bounded, blocking, multi-producer multi-consumer queue.
//
// Producers call push(), which blocks while the queue holds `capacity`
// items. Consumers call pop_wait(), which blocks until an item is available
// or the queue has been closed and drained. close() lets consumers finish
// without a separate completion signal.
template <typename T>
class ThreadSafeQueue
{
    std::queue<T> queue_;
    mutable std::mutex mutex_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    size_t capacity_;
    bool closed_ = false;

public:
    explicit ThreadSafeQueue(size_t capacity = std::numeric_limits<size_t>::max())
        : capacity_(capacity) {}
    ThreadSafeQueue(const ThreadSafeQueue<T> &) = delete;
    ThreadSafeQueue &operator=(const ThreadSafeQueue<T> &) = delete;
    ThreadSafeQueue(ThreadSafeQueue<T> &&) = delete;

    size_t size() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.size();
    }

    // Non-blocking pop. Returns an empty optional when the queue is empty.
    std::optional<T> pop()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (queue_.empty())
        {
            return std::nullopt;
        }
        std::optional<T> item(std::move(queue_.front()));
        queue_.pop();
        not_full_.notify_one();
        return item;
    }

    // Blocking pop. Returns an empty optional only once the queue has been
    // closed and every item has been consumed.
    std::optional<T> pop_wait()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        not_empty_.wait(lock, [this] { return !queue_.empty() || closed_; });
        if (queue_.empty())
        {
            return std::nullopt;
        }
        std::optional<T> item(std::move(queue_.front()));
        queue_.pop();
        not_full_.notify_one();
        return item;
    }

    // Blocking push. Waits while the queue is at capacity. Items pushed
    // after close() are dropped and false is returned.
    bool push(T &&item)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        not_full_.wait(lock, [this] { return queue_.size() < capacity_ || closed_; });
        if (closed_)
        {
            return false;
        }
        queue_.push(std::move(item));
        not_empty_.notify_one();
        return true;
    }

    bool push(const T &item)
    {
        T copy(item);
        return push(std::move(copy));
    }

    // Mark the queue as finished. Wakes all waiting consumers and producers.
    void close()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        closed_ = true;
        not_empty_.notify_all();
        not_full_.notify_all();
    }

    bool closed() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return closed_;
    }
};
