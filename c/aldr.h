/*
  Name:     aldr.h
  Purpose:  Fast sampling of random integers.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef ALDR_H
#define ALDR_H

// array
struct array_s
{
    int length;
    int *a;
};

// flattened ALDR tree
struct aldr_flat_s
{
    int length_breadths;
    int length_leaves_flat;
    int *breadths;
    int *leaves_flat;
};

void free_aldr_flat_s (struct aldr_flat_s x);
void free_array_s (struct array_s x);
struct aldr_flat_s preprocess_aldr_flat(int* a, int n);
struct aldr_flat_s preprocess_fldr_flat(int* a, int n);
struct array_s preprocess_aldr_enc(int* a, int n);
struct array_s preprocess_fldr_enc(int* a, int n);
int sample_aldr_flat(struct aldr_flat_s* f);
int sample_aldr_enc(struct array_s* x);
int bytes_sample_aldr_flat(struct aldr_flat_s *x);
int bytes_array(struct array_s *x);

#endif
