#include "hash_utils.h"
#include "configs.h"


int prev_a = 0;
int a_k = 0;
int a_k2 = 0;
int p = BIG_PRIME_NUMBER;

int zigma_hash(const char *seq_str, int length, int prev_hash, int seed, int gingle_gap) {
    int a = seed, res = 0;
    if (prev_hash == 0) {
        for (int i = 0; i < length + gingle_gap; i++) {
            if (i == length / 2)
                i += gingle_gap;
            res *= a;
            res += (seq_str[i] - 65);
            res %= p;
        }
        return res;
    } else {
        if (prev_a != a) {
            prev_a = a;
            a_k = a_k2 = 1;
            for (int i = 0; i < length - 1; i++)
                a_k = (a_k * a) % p;
            for (int i = 0; i < length / 2 - 1; i++)
                a_k2 = (a_k2 * a) % p;
        }
        res = prev_hash;
        res -= (seq_str[-1] - 65) * a_k;
        res += (seq_str[length / 2 - 1] - seq_str[length / 2 - 1 + gingle_gap]) * a_k2;
        res %= p;
        res += p;
        res %= p;
        res *= a;
        res += seq_str[length + gingle_gap - 1] - 65;
        res %= p;
        return res;
    }
}

