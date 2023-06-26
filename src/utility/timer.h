#pragma once

#include <chrono>

class Timer {
    typedef std::chrono::high_resolution_clock Clock_t;
    Clock_t::time_point _time;
public:
    inline void Reset() {
        _time = Clock_t::now();
    }
    inline float Elapsed() {
        auto time_now = Clock_t::now();
        float interval = (float)std::chrono::duration_cast<std::chrono::microseconds>(time_now - _time).count();
        _time = time_now;
        return interval / 1.0E6f;
    }
};
