#include "basket_min_hash.h"
#include "configs.h"

using namespace std;

BasketMinHash::BasketMinHash(const int sketch_window, const hash_function_t hash_function, const int max_hash_function)
        : sketch_window(sketch_window), hash_function(hash_function), max_hash_function(max_hash_function) {}

int *BasketMinHash::get_sketch(const Sequence &sequence, const unsigned int chunk_i, const int gingle_length,
                               const int gingle_gap) const {
    return get_sketch(sequence.seq_str, sequence.size, chunk_i, gingle_length, gingle_gap);
}

int *BasketMinHash::get_sketch(const char *seq_str, const unsigned long seq_str_len, const unsigned int chunk_i,
                               const int gingle_length, const int gingle_gap) const {
    unsigned int sketch_size = SKETCH_SIZE;
    auto sketch = new int[sketch_size];
    fill_n(sketch, sketch_size, max_hash_function);
    unsigned int sketch_step = max_hash_function / (sketch_size - 1);
    int new_hash[gingle_gap + 1] = {};
    for (int i = 0; i < seq_str_len - gingle_length - gingle_gap; i++)
        for (int j = 0; j <= gingle_gap; j++) {
            int new_hash_temp = new_hash[j] = hash_function(seq_str + i, gingle_length, new_hash[j], LOG_MAX_BASENUMBER,
                                                            j);
            unsigned int te = new_hash_temp / sketch_step;
            if (new_hash_temp < sketch[te])
                sketch[te] = new_hash_temp;
        }
    for (auto i = static_cast<int>(sketch_size - 2); i >= 0; i -= 1)
        if (sketch[i] == max_hash_function)
            sketch[i] = sketch[i + 1];//(i + 1) * sketch_step - 1;

    if (sketch_window > 1) {
        /*
            ssketch = [];
            for (int i =0; i <sketch_size; i+= sketch_window)
                ssketch.append(mmh3.hash(string(sketch + i:i + self.sketch_window])));
            return ssketch;
         */
    }
    return sketch;
}
