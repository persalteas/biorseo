#pragma once
#include <atomic>
#include <cassert>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>

class Pool
{

    private:
    std::queue<std::function<void()>> m_function_queue;
    std::mutex                        m_lock;
    std::condition_variable           m_data_condition;
    std::atomic<bool>                 m_accept_functions;

    public:
    Pool();
    ~Pool();
    std::mutex& get_mutex_reference(void);
    void push(std::function<void()> func);
    void done();
    void infinite_loop_func();
};

inline std::mutex &Pool::get_mutex_reference(void) { return m_lock; }

class quit_worker_exception : public std::exception {};