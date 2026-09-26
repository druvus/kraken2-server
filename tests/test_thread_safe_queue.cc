#include <doctest/doctest.h>

#include <atomic>
#include <thread>
#include <vector>

#include "thread_safe_queue.h"

TEST_CASE("ThreadSafeQueue basic push and pop")
{
    ThreadSafeQueue<int> q;
    CHECK(q.size() == 0);
    CHECK_FALSE(q.pop().has_value());
    q.push(1);
    q.push(2);
    CHECK(q.size() == 2);
    CHECK(q.pop().value() == 1);
    CHECK(q.pop_wait().value() == 2);
    CHECK_FALSE(q.pop().has_value());
}

TEST_CASE("ThreadSafeQueue moves items rather than copying")
{
    ThreadSafeQueue<std::unique_ptr<int>> q;  // move-only type compiles
    q.push(std::make_unique<int>(7));
    auto item = q.pop_wait();
    REQUIRE(item.has_value());
    CHECK(**item == 7);
}

TEST_CASE("ThreadSafeQueue pop_wait returns empty only after close and drain")
{
    ThreadSafeQueue<int> q;
    q.push(5);
    q.close();
    CHECK(q.closed());
    CHECK(q.pop_wait().value() == 5);
    CHECK_FALSE(q.pop_wait().has_value());
    // pushes after close are dropped
    CHECK_FALSE(q.push(6));
    CHECK_FALSE(q.pop_wait().has_value());
}

TEST_CASE("ThreadSafeQueue close wakes a blocked consumer")
{
    ThreadSafeQueue<int> q;
    std::atomic<bool> returned{false};
    std::thread consumer([&] {
        auto item = q.pop_wait();
        CHECK_FALSE(item.has_value());
        returned = true;
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    CHECK_FALSE(returned);
    q.close();
    consumer.join();
    CHECK(returned);
}

TEST_CASE("ThreadSafeQueue bounded push blocks until space is available")
{
    ThreadSafeQueue<int> q(2);
    q.push(1);
    q.push(2);
    std::atomic<bool> pushed{false};
    std::thread producer([&] {
        q.push(3);
        pushed = true;
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    CHECK_FALSE(pushed);
    CHECK(q.pop().value() == 1);
    producer.join();
    CHECK(pushed);
    CHECK(q.size() == 2);
}

TEST_CASE("ThreadSafeQueue delivers every item across producers and consumers")
{
    ThreadSafeQueue<int> q(8);
    const int per_producer = 500;
    std::atomic<long> sum{0};
    std::atomic<int> count{0};

    std::vector<std::thread> consumers;
    for (int c = 0; c < 3; c++)
        consumers.emplace_back([&] {
            while (auto item = q.pop_wait())
            {
                sum += *item;
                count++;
            }
        });
    std::vector<std::thread> producers;
    for (int p = 0; p < 4; p++)
        producers.emplace_back([&, p] {
            for (int i = 0; i < per_producer; i++)
                q.push(p * per_producer + i);
        });
    for (auto &t : producers) t.join();
    q.close();
    for (auto &t : consumers) t.join();

    long n = 4 * per_producer;
    CHECK(count == n);
    CHECK(sum == n * (n - 1) / 2);
}
