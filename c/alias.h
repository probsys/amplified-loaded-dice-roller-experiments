/*
  Name:     alias.h
  Purpose:  Exact weighted alias method, translated from Rust to C.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef ALIAS_H
#define ALIAS_H

#include <stdint.h>

// weighted alias index arrays
struct sample_weighted_alias_index_s {
    uint32_t length;
    uint32_t weight_sum;
    uint32_t *aliases;
    uint32_t *no_alias_odds;
};

void free_sample_weighted_alias_index(struct sample_weighted_alias_index_s x);
struct sample_weighted_alias_index_s preprocess_weighted_alias(int* a, int n);
uint32_t sample_weighted_alias_index(struct sample_weighted_alias_index_s *x);
int bytes_sample_weighted_alias_index(struct sample_weighted_alias_index_s *x);

#endif
