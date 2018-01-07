#include "hash_utils.h"
#include "configs.h"


int prev_a = 0;
int a_k = 0;
int a_k2 = 0;
int p = BIG_PRIME_NUMBER;
int length_2_1;

int zigma_hash(const char *seq_str, int length, int prev_hash, int shift_length, int gingle_gap) {
    int a = shift_length, res = 0;
    if (prev_hash == 0) {
        for (int i = 0; i < length + gingle_gap; i++) {
            if (i == length / 2)
                i += gingle_gap;
            res <<= a;
            res += (seq_str[i] - 65);
            res %= p;
        }
        return res;
    } else {
        if (prev_a != a) {
            length_2_1 = length / 2 - 1;
            prev_a = a;
            a_k = a_k2 = 1;
            for (int i = 0; i < length - 1; i++)
                a_k = (a_k << a) % p;
            for (int i = 0; i < length_2_1; i++)
                a_k2 = (a_k2 << a) % p;
        }
        res = prev_hash;
        res -= (seq_str[-1] - 65) * a_k;
        res += (seq_str[length_2_1] - seq_str[length_2_1 + gingle_gap]) * a_k2;
        res %= p;
        if (res < 0)
            res += p;
        res <<= a;
        res += seq_str[length + gingle_gap - 1] - 65;
        res %= p;
        return res;
    }
}

