#include <algorithm>
#include <memory>
#include "basket_min_hash.h"
#include "configs.h"

using namespace std;

BasketMinHash::BasketMinHash(const int sketch_window, const hash_function_t hash_function, const int max_hash_function)
        : sketch_window(sketch_window), hash_function(hash_function), max_hash_function(max_hash_function) {}

tuple<int *, char **>
BasketMinHash::get_sketch(const Sequence &sequence, const unsigned int chunk_i, const int gingle_length,
                          const int gingle_gap) const {
    return get_sketch(sequence.seq_str, sequence.size, chunk_i, gingle_length, gingle_gap);
}

//{'A': 904514536, 'C': 628428578, 'G': 631026863, 'T': 907208264, 'N': 161368370}
inline char nuc_acid_to_int(char nuc, char next_nuc='N') {
    if (next_nuc == 'G' || next_nuc == 'g') {
        if (nuc == 'T' || nuc == 't' || nuc == 'C' || nuc == 'c')
            return 5;  // TODO: Find a good value for this state
    }

    switch (nuc) {
        case 'T':
        case 't':
            return 3;
        case 'A':
        case 'a':
            return 2;
        case 'G':
        case 'g':
            return 1;
        case 'C':
        case 'c':
            return 0;
        case 'N':
        case 'n':
        case 'Y':
        case 'y':
        case 'U':
        case 'u':
        case 'K':
        case 'k':
        case 'W':
        case 'w':
        case 'B':
        case 'b':
        case 'D':
        case 'd':
        case 'H':
        case 'h':
            return 3;
        case 'R':
        case 'r':
        case 'M':
        case 'm':
        case 'V':
        case 'v':
            return 2;
        case 'S':
        case 's':
            return 1;
        default:
            return 3;
    }
}

tuple<int *, char **>
BasketMinHash::get_sketch(const char *seq_str, const unsigned long seq_str_len, const unsigned int chunk_i,
                          const int gingle_length, const int gingle_gap) const {
    unsigned int sketch_size = SKETCH_SIZE;
    unsigned int sketch_step = max_hash_function / (sketch_size - 1);

    auto sketch = new int[sketch_size];
    auto sketch_kmers = new char*[sketch_size];
    fill_n(sketch, sketch_size, max_hash_function);
    fill_n(sketch_kmers, sketch_size, NULL);

    // Change chars to int
    char seq_str_int[seq_str_len];
    for (int i = 0; i < seq_str_len; ++i) {
        char next_nuc = (i == seq_str_len - 1) ? 'N' : seq_str[i+1];
        seq_str_int[i] = nuc_acid_to_int(seq_str[i], next_nuc);
    }

    int new_hash[gingle_gap + 1] = {};
    for (int i = 0; i < seq_str_len - gingle_length - gingle_gap; i++) {
        for (int j = 0; j <= gingle_gap; j++) {
            int new_hash_temp = new_hash[j] = hash_function(seq_str_int + i, gingle_length, new_hash[j],
                                                            LOG_MAX_BASENUMBER, j);

            // Randomize hash
            new_hash_temp = (~new_hash_temp + (new_hash_temp << 21)) & BIG_PRIME_NUMBER;
            new_hash_temp = (new_hash_temp + (new_hash_temp << 3) + (new_hash_temp << 8)) & BIG_PRIME_NUMBER;
            new_hash_temp ^= (new_hash_temp >> 14);
            new_hash_temp = (new_hash_temp + (new_hash_temp << 2) + (new_hash_temp << 4)) & BIG_PRIME_NUMBER;

            // Find the basket this hash belongs to
            unsigned int sketch_i = new_hash_temp / sketch_step;
            // tmp->push_back(new_hash_temp);
            if (new_hash_temp < sketch[sketch_i])
                sketch[sketch_i] = new_hash_temp;

                if (sketch_kmers[sketch_i] == NULL)
                    sketch_kmers[sketch_i] = new char[gingle_length];
                std::copy(seq_str, seq_str + gingle_length, sketch_kmers[sketch_i]);
        }
    }

    // if there is no hash for a basket, fills it with the next basket
    for (auto i = static_cast<int>(sketch_size - 2); i >= 0; i -= 1) {
        if (sketch[i] == max_hash_function) {
            sketch[i] = sketch[i + 1];  // (i + 1) * sketch_step - 1;
            sketch_kmers[i] = sketch_kmers[i + 1];
        }
    }

    // sort(tmp->begin(), tmp->end());
    // if (tmp->size() < sketch_size)
    //    (*tmp)[0] = -1;
    // if (sketch_window > 1) {
    //     ssketch = [];
    //     for (int i =0; i <sketch_size; i+= sketch_window)
    //         ssketch.append(mmh3.hash(string(sketch + i:i + self.sketch_window])));
    //     return ssketch;
    // }
    // return &(*tmp)[0];

    return make_tuple(sketch, sketch_kmers);
}
