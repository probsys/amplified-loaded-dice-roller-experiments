/*
  Name:     flip.h
  Purpose:  Generating a sequence of pseudo-random bits.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef FLIP_H
#define FLIP_H

#include <stdint.h>

extern uint64_t NUM_RNG_CALLS;

extern uint32_t flip_k;
extern uint32_t flip_word;
extern uint32_t flip_pos;

uint32_t flip(void);
uint32_t flip_n(uint32_t n);
uint32_t uniform(uint32_t n);
uint32_t bernoulli(uint32_t numer, uint32_t denom);

// Macros for min and max.
#define max(a, b)            \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a, b)            \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

#endif
