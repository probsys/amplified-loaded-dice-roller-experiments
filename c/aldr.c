/*
  Name:     aldr.c
  Purpose:  Fast sampling of random integers.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>

#include "aldr.h"
#include "flip.h"

void free_aldr_flat_s (struct aldr_flat_s x) {
    free(x.breadths);
    free(x.leaves_flat);
}

void free_array_s (struct array_s x) {
    free(x.a);
}

struct aldr_flat_s preprocess_aldr_flat_k(int* a, int n, int kmul) {
    // assume k <= 31
    int m = 0;
    for (int i = 0; i < n; ++i) {
        m += a[i];
    }
    int k = 32 - __builtin_clz(m) - (1 == __builtin_popcount(m));
    int K = k * kmul;                       // depth
    long long c = (1ll << K) / m;           // amplification factor
    long long r = (1ll << K) % m;           // reject weight

    int num_leaves = __builtin_popcountll(r);
    for (int i = 0; i < n; ++i) {
        num_leaves += __builtin_popcountll(c * a[i]);
    }

    int *breadths = calloc(K + 1, sizeof(int));
    int *leaves_flat = calloc(num_leaves, sizeof(int));

    int location = 0;
    for(int j = 0; j <= K; j++) {
        long long bit = (1ll << (K - j));
        if (r & bit) {
            leaves_flat[location] = 0;
            ++breadths[j];
            ++location;
        }
        for (int i = 0; i < n; ++i) {
            long long Qi = c*a[i];
            if (Qi & bit) {
                leaves_flat[location] = i + 1;
                ++breadths[j];
                ++location;
            }
        }
    }

    return (struct aldr_flat_s){
            .length_breadths = K+1,
            .length_leaves_flat = num_leaves,
            .breadths = breadths,
            .leaves_flat = leaves_flat
        };
}

struct aldr_flat_s preprocess_aldr_flat(int* a, int n) {
    return preprocess_aldr_flat_k(a, n, 2);
}

struct aldr_flat_s preprocess_fldr_flat(int* a, int n) {
    return preprocess_aldr_flat_k(a, n, 1);
}

struct array_s preprocess_aldr_enc_k(int* a, int n, int kmul) {
    // assume k <= 31 like fldr does
    int m = 0;
    for (int i = 0; i < n; ++i) {
        m += a[i];
    }
    int k = 32 - __builtin_clz(m) - (1 == __builtin_popcount(m));
    int K = k * kmul;
    long long c = (1ll << K) / m;
    long long r = (1ll << K) % m;

    // flattened but 50% sparse encoding
    int num_leaves = __builtin_popcountll(r);
    for (int i = 0; i < n; ++i) {
        num_leaves += __builtin_popcountll(c * a[i]);
    }
    int *enc = calloc((num_leaves << 1) - 1, sizeof(int));
    int prev_location = 0;
    int prev_length = 1;
    int location = 0;
    for(int j = 0; j <= K; j++) {
        long long bit = (1ll << (K - j));
        int next_length = prev_length;
        if (r & bit) {
            // reject outcome: flip and go to child of root
            enc[location++] = 1;
        }
        if (r & bit) {
            enc[location++] = ~0;
        }
        for (int i = 0; i < n; ++i) {
            if ((c*a[i]) & bit) {
                enc[location++] = ~(i+1);
            }
        }
        for ( ; location < prev_length; ++location) {
            enc[location] = next_length;
            next_length += 2;
        }
        prev_length = next_length;
    }

    return (struct array_s){ .length = (num_leaves << 1) - 1, .a = enc };
}

struct array_s preprocess_aldr_enc(int* a, int n) {
    return preprocess_aldr_enc_k(a, n, 2);
}

struct array_s preprocess_fldr_enc(int* a, int n) {
    return preprocess_aldr_enc_k(a, n, 1);
}

int sample_aldr_flat(struct aldr_flat_s* f) {
    while (1) {
        int depth = 0;
        int location = 0;
        int val = 0;
        for (;;) {
            if (val < f->breadths[depth]) {
                int ans = f->leaves_flat[location + val];
                if (ans) return ans - 1;
                else break;
            }
            location += f->breadths[depth];
            val = ((val - f->breadths[depth]) << 1) | flip();
            ++depth;
        }
    }
}

int sample_aldr_enc(struct array_s* x) {
    int c = x->a[0];
    for (;;) {
        if (c < 0) {
            // note that this implementation uses all pointers
            // instead of "n" with special reject logic
            return (~c) - 1;
        }
        c = x->a[c+flip()];
    }
}

int bytes_sample_aldr_flat(struct aldr_flat_s *x) {
    // this doesn't count the length variables themselves
    // because we don't need them and just added them here
    // for easy byte counting
    return
        x->length_breadths * sizeof(x->breadths[0])
            + x->length_leaves_flat * sizeof(x->leaves_flat[0]);
}

int bytes_array(struct array_s *x) {
    return x->length * sizeof(x->a[0]) + sizeof(x->length);
}
