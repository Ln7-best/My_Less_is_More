#ifndef READERWRITERQUEUE_H
#define READERWRITERQUEUE_H
#include <cstddef>
#include <cstdio>
#include <vector>
#pragma once

#define NDEBUG
#include "atomicops.h"
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <emmintrin.h>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

#ifndef EXCEPTIONS_ENABLED
#if (defined(_MSC_VER) && defined(_CPPUNWIND)) ||                              \
    (defined(__GNUC__) && defined(__EXCEPTIONS)) ||                            \
    (!defined(_MSC_VER) && !defined(__GNUC__))
#define EXCEPTIONS_ENABLED
#endif
#endif

#ifndef HAS_EMPLACE
#if !defined(_MSC_VER) ||                                                      \
    _MSC_VER >=                                                                \
        1800 // variadic templates: either a non-MS compiler or VS >= 2013
#define HAS_EMPLACE 1
#endif
#endif


#ifndef MAYBE_ALIGN_TO_CACHELINE
#define MAYBE_ALIGN_TO_CACHELINE ALIGN(CACHE_LINE_SIZE)
#endif

template <typename T, size_t MAX_BLOCK_SIZE = 64>
class MAYBE_ALIGN_TO_CACHELINE ReaderWriterQueue {

public:
  typedef T value_type;
  size_t expansion_cnt;
  static size_t queue_size;
  std::vector<double> make_block_time;

  NO_TSAN explicit ReaderWriterQueue(size_t size = 32)
#ifndef NDEBUG
      : enqueuing(false), dequeuing(false)
#endif
  {
    expansion_cnt = 0;
    assert(MAX_BLOCK_SIZE == ceilToPow2(MAX_BLOCK_SIZE) &&
           "MAX_BLOCK_SIZE must be a power of 2");
    assert(MAX_BLOCK_SIZE >= 2 && "MAX_BLOCK_SIZE must be at least 2");

    Block *firstBlock = nullptr;

    largestBlockSize = ceilToPow2(
        size + 1); 
    if (largestBlockSize > MAX_BLOCK_SIZE * 2) {
      size_t initialBlockCount =
          (size + MAX_BLOCK_SIZE * 2 - 3) / (MAX_BLOCK_SIZE - 1);
      largestBlockSize = MAX_BLOCK_SIZE;
      Block *lastBlock = nullptr;
      for (size_t i = 0; i != initialBlockCount; ++i) {
        auto block = make_block(largestBlockSize);
        if (block == nullptr) {
#ifdef EXCEPTIONS_ENABLED
          throw std::bad_alloc();
#else
          abort();
#endif
        }
        if (firstBlock == nullptr) {
          firstBlock = block;
        } else {
          lastBlock->next = block;
        }
        lastBlock = block;
        block->next = firstBlock;
      }
    } else {
      firstBlock = make_block(largestBlockSize);
      if (firstBlock == nullptr) {
#ifdef EXCEPTIONS_ENABLED
        throw std::bad_alloc();
#else
        abort();
#endif
      }
      firstBlock->next = firstBlock;
    }
    frontBlock = firstBlock;
    tailBlock = firstBlock;
    fence(memory_order_sync);
  }

  NO_TSAN ReaderWriterQueue(ReaderWriterQueue &&other)
      : frontBlock(other.frontBlock.load()), tailBlock(other.tailBlock.load()),
        largestBlockSize(other.largestBlockSize)
#ifndef NDEBUG
        ,
        enqueuing(false), dequeuing(false)
#endif
  {
    other.largestBlockSize = 32;
    Block *b = other.make_block(other.largestBlockSize);
    if (b == nullptr) {
#ifdef EXCEPTIONS_ENABLED
      throw std::bad_alloc();
#else
      abort();
#endif
    }
    b->next = b;
    other.frontBlock = b;
    other.tailBlock = b;
  }

  ReaderWriterQueue &operator=(ReaderWriterQueue &&other) NO_TSAN {
    Block *b = frontBlock.load();
    frontBlock = other.frontBlock.load();
    other.frontBlock = b;
    b = tailBlock.load();
    tailBlock = other.tailBlock.load();
    other.tailBlock = b;
    std::swap(largestBlockSize, other.largestBlockSize);
    return *this;
  }

  NO_TSAN ~ReaderWriterQueue() {
    fence(memory_order_sync);

    Block *frontBlock_ = frontBlock;
    Block *block = frontBlock_;
    do {
      Block *nextBlock = block->next;
      size_t blockFront = block->front;
      size_t blockTail = block->tail;

      for (size_t i = blockFront; i != blockTail;
           i = (i + 1) & block->sizeMask) {
        auto element = reinterpret_cast<T *>(block->data + i * sizeof(T));
        element->~T();
        (void)element;
      }

      auto rawBlock = block->rawThis;
      block->~Block();
      std::free(rawBlock);
      block = nextBlock;
    } while (block != frontBlock_);
  }

  inline bool enqueue(T const &element) NO_TSAN {
    return inner_enqueue(element);
  }

  inline bool enqueue(T &&element) NO_TSAN {
    return inner_enqueue(std::forward<T>(element));
  }

  inline bool enqueue_fast(T const &element) NO_TSAN {
    return inner_enqueue_fast(element);
  }

  inline bool enqueue_fast(T &&element) NO_TSAN {
    return inner_enqueue_fast(std::forward<T>(element));
  }
#if HAS_EMPLACE
#endif

  template <typename U> bool try_dequeue(U &result) NO_TSAN {
#ifndef NDEBUG
    ReentrantGuard guard(this->dequeuing);
#endif

    Block *frontBlock_ = frontBlock.load();
    size_t blockTail = frontBlock_->localTail;
    size_t blockFront = frontBlock_->front.load();

    if (blockFront != blockTail ||
        blockFront != (frontBlock_->localTail = frontBlock_->tail.load())) {
      fence(memory_order_acquire);

    non_empty_front_block:
      auto element =
          reinterpret_cast<T *>(frontBlock_->data + blockFront * sizeof(T));
      result = std::move(*element);
      blockFront = (blockFront + 1) & frontBlock_->sizeMask;

      fence(memory_order_release);
      frontBlock_->front = blockFront;
    } else if (frontBlock_ != tailBlock.load()) {
      fence(memory_order_acquire);

      frontBlock_ = frontBlock.load();
      blockTail = frontBlock_->localTail = frontBlock_->tail.load();
      blockFront = frontBlock_->front.load();
      fence(memory_order_acquire);

      if (blockFront != blockTail) {
        goto non_empty_front_block;
      }

      Block *nextBlock = frontBlock_->next;

      size_t nextBlockFront = nextBlock->front.load();
      size_t nextBlockTail = nextBlock->localTail = nextBlock->tail.load();
      fence(memory_order_acquire);


      fence(memory_order_release); 
      frontBlock = frontBlock_ = nextBlock;

      compiler_fence(memory_order_release); // Not strictly needed

      auto element =
          reinterpret_cast<T *>(frontBlock_->data + nextBlockFront * sizeof(T));

      result = std::move(*element);
      nextBlockFront = (nextBlockFront + 1) & frontBlock_->sizeMask;

      fence(memory_order_release);
      frontBlock_->front = nextBlockFront;
    } else {
      return false;
    }

    return true;
  }

  inline size_t max_capacity() const {
    size_t result = 0;
    Block *frontBlock_ = frontBlock.load();
    Block *block = frontBlock_;
    do {
      fence(memory_order_acquire);
      result += block->sizeMask;
      block = block->next.load();
    } while (block != frontBlock_);
    return result;
  }

#if HAS_EMPLACE
  template <typename... Args>
  bool inner_enqueue(Args &&...args) NO_TSAN
#else
  template <typename U>
  bool inner_enqueue(U &&element) NO_TSAN
#endif
  {
#ifndef NDEBUG
    ReentrantGuard guard(this->enqueuing);
#endif


    Block *tailBlock_ = tailBlock.load();
    size_t blockFront = tailBlock_->localFront;
    size_t blockTail = tailBlock_->tail.load();

    size_t nextBlockTail = (blockTail + 1) & tailBlock_->sizeMask;
    if (nextBlockTail != blockFront ||
        nextBlockTail != (tailBlock_->localFront = tailBlock_->front.load())) {
      fence(memory_order_acquire);
      char *location = tailBlock_->data + blockTail * sizeof(T);
#if HAS_EMPLACE
      new (location) T(std::forward<Args>(args)...);
#else
      new (location) T(std::forward<U>(element));
#endif

      fence(memory_order_release);
      tailBlock_->tail = nextBlockTail;
    } else {
      fence(memory_order_acquire);
      if (tailBlock_->next.load() != frontBlock) {

        fence(memory_order_acquire); 
        Block *tailBlockNext = tailBlock_->next.load();
        size_t nextBlockFront = tailBlockNext->localFront =
            tailBlockNext->front.load();
        nextBlockTail = tailBlockNext->tail.load();
        fence(memory_order_acquire);

        assert(nextBlockFront == nextBlockTail);
        tailBlockNext->localFront = nextBlockFront;

        char *location = tailBlockNext->data + nextBlockTail * sizeof(T);
#if HAS_EMPLACE
        new (location) T(std::forward<Args>(args)...);
#else
        new (location) T(std::forward<U>(element));
#endif

        tailBlockNext->tail = (nextBlockTail + 1) & tailBlockNext->sizeMask;

        fence(memory_order_release);
        tailBlock = tailBlockNext;
      } else {

        auto newBlockSize = largestBlockSize >= MAX_BLOCK_SIZE
                                ? largestBlockSize
                                : largestBlockSize * 2;
        auto start_time = std::chrono::high_resolution_clock::now();
        auto newBlock = make_block(newBlockSize);
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration =
            end_time - start_time;
        make_block_time.emplace_back(duration.count());
        if (newBlock == nullptr) {
          return false;
        }
        expansion_cnt++;
        largestBlockSize = newBlockSize;

#if HAS_EMPLACE
        new (newBlock->data) T(std::forward<Args>(args)...);
#else
        new (newBlock->data) T(std::forward<U>(element));
#endif
        assert(newBlock->front == 0);
        newBlock->tail = newBlock->localTail = 1;

        newBlock->next = tailBlock_->next.load();
        tailBlock_->next = newBlock;

        fence(memory_order_release);
        tailBlock = newBlock;
        return false;
      }
    }

    return true;
  }
  #if HAS_EMPLACE
  template <typename... Args>
  bool inner_enqueue_fast(Args &&...args) NO_TSAN
#else
  template <typename U>
  bool inner_enqueue_fast(U &&element) NO_TSAN
#endif
  {
#ifndef NDEBUG
    ReentrantGuard guard(this->enqueuing);
#endif
    Block *tailBlock_ = tailBlock.load();
    size_t blockFront = tailBlock_->localFront;
    size_t blockTail = tailBlock_->tail.load();
    size_t nextBlockTail = (blockTail + 1) & tailBlock_->sizeMask;
    if (nextBlockTail != blockFront ||
        nextBlockTail != (tailBlock_->localFront = tailBlock_->front.load())) {
      fence(memory_order_acquire);
      // This block has room for at least one more element
      char *location = tailBlock_->data + blockTail * sizeof(T);
#if HAS_EMPLACE
      new (location) T(std::forward<Args>(args)...);
#else
      new (location) T(std::forward<U>(element));
#endif

      fence(memory_order_release);
      tailBlock_->tail = nextBlockTail;
    } else {
      fence(memory_order_acquire);
      if (tailBlock_->next.load() != frontBlock) {

        fence(memory_order_acquire); 

        Block *tailBlockNext = tailBlock_->next.load();
        size_t nextBlockFront = tailBlockNext->localFront =
            tailBlockNext->front.load();
        nextBlockTail = tailBlockNext->tail.load();
        fence(memory_order_acquire);

        assert(nextBlockFront == nextBlockTail);
        tailBlockNext->localFront = nextBlockFront;

        char *location = tailBlockNext->data + nextBlockTail * sizeof(T);
#if HAS_EMPLACE
        new (location) T(std::forward<Args>(args)...);
#else
        new (location) T(std::forward<U>(element));
#endif

        tailBlockNext->tail = (nextBlockTail + 1) & tailBlockNext->sizeMask;

        fence(memory_order_release);
        tailBlock = tailBlockNext;
      } else {
        return false;
      }
    }

    return true;
  }


  ReaderWriterQueue(ReaderWriterQueue const &) {}

  ReaderWriterQueue &operator=(ReaderWriterQueue const &) {}

  FORCEINLINE static size_t ceilToPow2(size_t x) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    for (size_t i = 1; i < sizeof(size_t); i <<= 1) {
      x |= x >> (i << 3);
    }
    ++x;
    return x;
  }

  template <typename U>
  static FORCEINLINE char *align_for(char *ptr) NO_TSAN {
    const std::size_t alignment = std::alignment_of<U>::value;
    return ptr +
           (alignment - (reinterpret_cast<std::uintptr_t>(ptr) % alignment)) %
               alignment;
  }

  struct Block {

    weak_atomic<size_t> front; 
    size_t
        localTail; 

    char cachelineFiller0[CACHE_LINE_SIZE -
                          sizeof(weak_atomic<size_t>) - sizeof(size_t)];
    weak_atomic<size_t> tail; 
    size_t localFront;

    char cachelineFiller1[CACHE_LINE_SIZE -
                          sizeof(weak_atomic<size_t>) -
                          sizeof(size_t)]; 
    weak_atomic<Block *> next;            

    char *data; 

    const size_t sizeMask;

    NO_TSAN Block(size_t const &_size, char *_rawThis, char *_data)
        : front(0UL), localTail(0), tail(0UL), localFront(0), next(nullptr),
          data(_data), sizeMask(_size - 1), rawThis(_rawThis) {}

  private:
    Block &operator=(Block const &);

  public:
    char *rawThis;
  };

  static Block *make_block(size_t capacity) NO_TSAN {
    auto size = sizeof(Block) + std::alignment_of<Block>::value - 1;
    size += sizeof(T) * capacity + std::alignment_of<T>::value - 1;
    auto newBlockRaw = static_cast<char *>(std::malloc(size));
    volatile char *touch = newBlockRaw;
    for (size_t i = 0; i < size; i += 4096) {
        touch[i] = 0; 
    }
    queue_size += size;
    if (newBlockRaw == nullptr) {
      return nullptr;
    }
    auto newBlockAligned = align_for<Block>(newBlockRaw);
    auto newBlockData = align_for<T>(newBlockAligned + sizeof(Block));
    return new (newBlockAligned) Block(capacity, newBlockRaw, newBlockData);
  }

  weak_atomic<Block *>
      frontBlock; 
  std::atomic<int64_t> counts{0};

  char cachelineFiller[CACHE_LINE_SIZE -
                       sizeof(weak_atomic<Block *>)];
  weak_atomic<Block *>
      tailBlock;

  size_t largestBlockSize;
#ifndef NDEBUG
  weak_atomic<bool> enqueuing;
  mutable weak_atomic<bool> dequeuing;
#endif
};
template <typename T, size_t MAX_BLOCK_SIZE>
size_t ReaderWriterQueue<T, MAX_BLOCK_SIZE>::queue_size = 0;


#endif // PTHREAD_READERWRITERQUEUE_H
