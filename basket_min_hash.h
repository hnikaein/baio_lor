#ifndef C_LONG_READ_ALIGNER_BASKET_MIN_HASH_H
#define C_LONG_READ_ALIGNER_BASKET_MIN_HASH_H

#include <functional>
#include <vector>
#include "bio_utils/sequence.h"


typedef std::function<int(const char *const, const int, const int, const int, const int)> hash_function_t;


class BasketMinHash {
public:
    BasketMinHash(int sketch_ratio, int sketch_window, hash_function_t hash_function, int max_hash_function);

    int *get_sketch(const Sequence &sequence, unsigned int chunk_size, int gingle_length, int gingle_gap) const;

    int *get_sketch(const char *seq_str, unsigned long seq_str_len, unsigned int chunk_size, int gingle_length,
                    int gingle_gap) const;

    int sketch_ratio;
private:
    int sketch_window, max_hash_function;
    hash_function_t hash_function;
};

#endif //C_LONG_READ_ALIGNER_BASKET_MIN_HASH_H
