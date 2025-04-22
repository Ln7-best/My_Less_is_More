#ifndef ATOMICOPS_H
#define ATOMICOPS_H
#pragma once

#include <cerrno>
#include <cassert>
#include <type_traits>
#include <cerrno>
#include <cstdint>
#include <ctime>



#ifndef NO_TSAN
#define NO_TSAN
#endif
#ifndef TSAN_ANNOTATE_RELEASE
#define TSAN_ANNOTATE_RELEASE()
#define TSAN_ANNOTATE_ACQUIRE()
#endif



#define FORCEINLINE inline
#define ALIGN(x) __attribute__((aligned(x)))


#include <atomic>
#include <utility>

enum memory_order {
    memory_order_relaxed,
    memory_order_acquire,
    memory_order_release,
    memory_order_acq_rel,
    memory_order_seq_cst,
    memory_order_sync = memory_order_seq_cst
};

inline void compiler_fence(memory_order order) NO_TSAN
{
    switch (order) {
        case memory_order_relaxed: break;
        case memory_order_acquire: std::atomic_signal_fence(std::memory_order_acquire); break;
        case memory_order_release: std::atomic_signal_fence(std::memory_order_release); break;
        case memory_order_acq_rel: std::atomic_signal_fence(std::memory_order_acq_rel); break;
        case memory_order_seq_cst: std::atomic_signal_fence(std::memory_order_seq_cst); break;
        default: assert(false);
    }
}

inline void fence(memory_order order) NO_TSAN
{
    switch (order) {
        case memory_order_relaxed: break;
        case memory_order_acquire: TSAN_ANNOTATE_ACQUIRE(); std::atomic_thread_fence(std::memory_order_acquire); break;
        case memory_order_release: TSAN_ANNOTATE_RELEASE(); std::atomic_thread_fence(std::memory_order_release); break;
        case memory_order_acq_rel: TSAN_ANNOTATE_ACQUIRE(); TSAN_ANNOTATE_RELEASE(); std::atomic_thread_fence(std::memory_order_acq_rel); break;
        case memory_order_seq_cst: TSAN_ANNOTATE_ACQUIRE(); TSAN_ANNOTATE_RELEASE(); std::atomic_thread_fence(std::memory_order_seq_cst); break;
        default: assert(false);
    }
}

template<typename T>
class weak_atomic
{
public:
    NO_TSAN weak_atomic() : value() { }

    template<typename U> NO_TSAN weak_atomic(U&& x) : value(std::forward<U>(x)) {  }

    NO_TSAN weak_atomic(weak_atomic const& other) : value(other.load()) {  }
    NO_TSAN weak_atomic(weak_atomic&& other) : value(std::move(other.load())) {  }

    inline operator T() const NO_TSAN { return load(); }

    template<typename U>
    inline weak_atomic const& operator=(U&& x) NO_TSAN
    {
        value.store(std::forward<U>(x), std::memory_order_relaxed);
        return *this;
    }

    inline weak_atomic const& operator=(weak_atomic const& other) NO_TSAN
    {
        value.store(other.value.load(std::memory_order_relaxed), std::memory_order_relaxed);
        return *this;
    }

    inline T load() const NO_TSAN { return value.load(std::memory_order_relaxed); }

private:
    std::atomic<T> value;
};

#endif

