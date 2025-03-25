/*
  Name:     alias.c
  Purpose:  Exact weighted alias method, translated from Rust to C.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include "flip.h"
#include "alias.h"

void free_sample_weighted_alias_index(struct sample_weighted_alias_index_s x) {
    free(x.aliases);
    free(x.no_alias_odds);
}

struct Aliases {
    uint32_t *aliases;
    uint32_t smalls_head;
    uint32_t bigs_head;
};

/// This struct is designed to contain three data structures at once,
/// sharing the same memory. More precisely it contains two linked lists
/// and an alias map, which will be the output of this method. To keep
/// the three data structures from getting in each other's way, it must
/// be ensured that a single index is only ever in one of them at the
/// same time.
struct Aliases aliases_new(uint32_t n) {
    return (struct Aliases) {
        .aliases = calloc(n, sizeof(uint32_t)),
        .smalls_head = UINT32_MAX,
        .bigs_head = UINT32_MAX
    };
}

void push_small(struct Aliases *aliases, uint32_t idx) {
    aliases->aliases[idx] = aliases->smalls_head;
    aliases->smalls_head = idx;
}

void push_big(struct Aliases *aliases, uint32_t idx) {
    aliases->aliases[idx] = aliases->bigs_head;
    aliases->bigs_head = idx;
}

uint32_t pop_small(struct Aliases *aliases) {
    uint32_t idx = aliases->smalls_head;
    aliases->smalls_head = aliases->aliases[idx];
    return idx;
}

uint32_t pop_big(struct Aliases *aliases) {
    uint32_t idx = aliases->bigs_head;
    aliases->bigs_head = aliases->aliases[idx];
    return idx;
}

bool smalls_is_empty(struct Aliases *aliases) {
    return aliases->smalls_head == UINT32_MAX;
}

bool bigs_is_empty(struct Aliases *aliases) {
    return aliases->bigs_head == UINT32_MAX;
}

void set_alias(struct Aliases *aliases, uint32_t idx, uint32_t alias) {
    aliases->aliases[idx] = alias;
}

/// Creates a new [`WeightedAliasIndex`].
///
/// Returns an error if:
/// - The vector is empty.
/// - The vector is longer than `u32::MAX`.
/// - For any weight `w`: `w < 0` or `w > max` where `max = W::MAX /
///   weights.len()`.
/// - The sum of weights is zero.
struct sample_weighted_alias_index_s preprocess_weighted_alias(int* a, int n) {
    assert(n > 0);
    assert(n < UINT32_MAX);
    uint32_t max_weight_size = UINT32_MAX / n;
    for (uint32_t i = 0; i < n; ++i) {
        assert(0 <= a[i]);
        assert(a[i] <= max_weight_size);
    }

    // The sum of weights will represent 100% of no alias odds.
    uint32_t weight_sum = 0;
    for (uint32_t i = 0; i < n; ++i) {
        weight_sum += a[i];
    }
    assert(weight_sum >= 0);

    uint32_t *no_alias_odds = calloc(n, sizeof(uint32_t));
    for (uint32_t i = 0; i < n; ++i) {
        no_alias_odds[i] = a[i] * n;
    }

    struct Aliases aliases = aliases_new(n);

    // Split indices into those with small weights and those with big weights.
    for (uint32_t i = 0; i < n; ++i) {
        if (no_alias_odds[i] < weight_sum) {
            push_small(&aliases, i);
        } else {
            push_big(&aliases, i);
        }
    }

    // Build the alias map by finding an alias with big weight for each index with
    // small weight.
    while (!smalls_is_empty(&aliases) && !bigs_is_empty(&aliases)) {
        uint32_t small = pop_small(&aliases);
        uint32_t big = pop_big(&aliases);
        set_alias(&aliases, small, big);
        no_alias_odds[big] -= weight_sum - no_alias_odds[small];
        if (no_alias_odds[big] < weight_sum) {
            push_small(&aliases, big);
        } else {
            push_big(&aliases, big);
        }
    }

    // The remaining indices should have no alias odds of about 100%. This is due to
    // numeric accuracy. Otherwise they would be exactly 100%.
    while (!smalls_is_empty(&aliases)) {
        no_alias_odds[pop_small(&aliases)] = weight_sum;
    }
    while (!bigs_is_empty(&aliases)) {
        no_alias_odds[pop_big(&aliases)] = weight_sum;
    }

    return (struct sample_weighted_alias_index_s) {
        .length = n,
        .weight_sum = weight_sum,
        .aliases = aliases.aliases,
        .no_alias_odds = no_alias_odds
    };
}

uint32_t sample_weighted_alias_index(struct sample_weighted_alias_index_s *x) {
    uint32_t uniform_index = uniform(x->length);
    if (bernoulli(x->no_alias_odds[uniform_index], x->weight_sum)) {
        return uniform_index;
    } else {
        return x->aliases[uniform_index];
    }
}

int bytes_sample_weighted_alias_index(struct sample_weighted_alias_index_s *x) {
    return
        x->length * sizeof(x->aliases[0])
            + x->length * sizeof(x->no_alias_odds[0])
            + sizeof(x->length)
            + sizeof(x->weight_sum);
}
