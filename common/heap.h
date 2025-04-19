#ifndef HEAP_H
#define HEAP_H
#define MINUPDATE_THRESHOLD 10
#include <algorithm>
#include <unordered_map>
#include <atomic>
#include "cuckooMap.h"
template<typename DATA_TYPE, typename COUNT_TYPE>
class Heap{
public:
    std::atomic<int32_t> min_val;
    // int32_t min_val;
    struct KV{
        DATA_TYPE key;
        COUNT_TYPE value;
    };
    uint64_t update_round;
    typedef CuckooMap<DATA_TYPE, uint32_t> Cuckoo;
    typedef std::unordered_map<DATA_TYPE, COUNT_TYPE> HashMap;
    typedef typename Cuckoo::Entry Entry;
    std::atomic<bool> full_flag;
    bool first_update;
    Heap(uint32_t _SIZE){
        SIZE = _SIZE;
        update_round = 0;
        first_update = false;
        full_flag.store(false);
        mp = new Cuckoo(SIZE);
        heap = new KV[SIZE];
        memset(heap, 0, sizeof(KV) * SIZE);
        min_val = 0;
    }

    ~Heap(){
        delete mp;
        delete [] heap;
    }

    static uint32_t Size2Memory(uint32_t size){
        return size * ((sizeof(DATA_TYPE) + sizeof(uint32_t)) / LOAD
                       + sizeof(DATA_TYPE) + sizeof(COUNT_TYPE));
    }

    static uint32_t Memory2Size(uint32_t memory){
        return memory / ((sizeof(DATA_TYPE) + sizeof(uint32_t)) / LOAD
                         + sizeof(DATA_TYPE) + sizeof(COUNT_TYPE));
    }

    inline COUNT_TYPE min(){
    if(full_flag.load(std::memory_order_relaxed))
        return min_val.load(std::memory_order_relaxed);
    else
        return 0;
        
    }

    void Insert(const DATA_TYPE item, const COUNT_TYPE frequency){
        if(this->isFull()){
            uint64_t origin_min = heap[0].value;
            if(frequency > heap[0].value){
                Entry temp = mp->Lookup(item);
                if(temp.bucket != nullptr){
                    uint32_t pos = temp.bucket->values[temp.position];
                    if(heap[pos].value<frequency)
                        heap[pos].value = frequency;
                    Heap_Down(temp, pos);
                }
                else{
                    mp->Delete(heap[0].key);
                    heap[0].value = frequency;
                    heap[0].key = item;
                    mp->Insert(item, 0);
                    Entry temp = mp->Lookup(item);
                    this->Heap_Down(temp, 0);
                }
                // if(update_round%MINUPDATE_THRESHOLD==0)
                //     min_val.store(heap[0].value);
                if(origin_min != heap[0].value&&update_round%MINUPDATE_THRESHOLD==0)
                {
                    // if(update_round % MINUPDATE_THRESHOLD)
                // if(origin_min != heap[0].value&&update_round%MINUPDATE_THRESHOLD==0)
                    min_val.store(heap[0].value);
                }
                update_round++;

            }
            return;
        }

        Entry temp = mp->Lookup(item);
        if(temp.bucket != nullptr){
            uint32_t pos = temp.bucket->values[temp.position];
            if(heap[pos].value<frequency)
                heap[pos].value = frequency;
            Heap_Down(temp, pos);
        }
        else{
            uint32_t pos = mp->size();
            heap[pos].value = frequency;
            heap[pos].key = item;
            mp->Insert(item, pos);
            Heap_Up(pos);
        }
        if(!first_update && this->isFull())
        {
            full_flag.store(true);
            first_update = true;
            min_val.store(heap[0].value);
        }
    }

    //COUNT_TYPE Query(const DATA_TYPE& item){
    //    return mp->Lookup(item)?  heap[(*mp)[item]].value : 0;
    //}

    HashMap AllQuery(){
        HashMap ret;
        uint32_t size = mp->size();
        for(uint32_t i = 0;i < size;++i){
            ret[heap[i].key] = heap[i].value;
        }
        return ret;
    }

//protected:
    uint32_t SIZE;
    Cuckoo* mp;
    KV* heap;

    inline bool isFull(){
        return mp->size() >= SIZE;
    }

    void Heap_Down(Entry temp, uint32_t pos){
        uint32_t upper = mp->size();

        while (pos < upper / 2) {
            uint32_t left = 2 * pos + 1, right = 2 * pos + 2;
            uint32_t replace = pos;

            if (left < upper && heap[left].value < heap[replace].value)
                replace = left;
            if (right < upper && heap[right].value< heap[replace].value)
                replace = right;

            if (replace != pos) {
                temp.bucket->values[temp.position] = replace;
                KV temp = heap[pos];
                heap[pos] = heap[replace];
                heap[replace] = temp;
                mp->Replace(heap[pos].key, pos);
                pos = replace;
            }
            else return;
        }
    }

    void Heap_Up(uint32_t pos) {
        while (pos > 1) {
            uint32_t parent = (pos - 1) / 2;
            if (heap[parent].value <= heap[pos].value)
                break;

            KV temp = heap[pos];
            heap[pos] = heap[parent];
            heap[parent] = temp;
            mp->Replace(heap[pos].key, pos);
            mp->Replace(heap[parent].key, parent);
            pos = parent;
        }
    }
};

#endif
