#include <doctest/doctest.h>

#include <atomic>
#include <chrono>
#include <set>
#include <mutex>

#include "worker_pool.h"

TEST_CASE("WorkerPool runs every queued task before destruction")
{
    std::atomic<int> count{0};
    {
        WorkerPool pool(3);
        CHECK(pool.thread_count() == 3);
        for (int i = 0; i < 200; i++)
            pool.push([&] { count++; });
    }
    CHECK(count == 200);
}

TEST_CASE("WorkerPool with zero threads requested uses hardware concurrency")
{
    WorkerPool pool(0);
    CHECK(pool.thread_count() >= 1);
}

TEST_CASE("WorkerPool runs tasks on more than one thread")
{
    std::mutex m;
    std::set<std::thread::id> ids;
    {
        WorkerPool pool(4);
        for (int i = 0; i < 16; i++)
            pool.push([&] {
                std::this_thread::sleep_for(std::chrono::milliseconds(5));
                std::lock_guard<std::mutex> lock(m);
                ids.insert(std::this_thread::get_id());
            });
    }
    CHECK(ids.size() > 1);
}
