/*
  Name:     main.c
  Purpose:  Main function for benchmarking ALDR and the alias method.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "flip.h"
#include "macros.c"
#include "aldr.h"
#include "alias.h"

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("usage: %s sampler path\n", argv[0]);
        exit(0);
    }
    char *sampler = argv[1];
    char *path = argv[2];

    // check if sampler ends in ".osrng".
    // If so, set flip_k = 32 and remove ".osrng" from sampler.
    if (strstr(sampler, ".osrng") != NULL) {
        flip_k = 32;
        sampler[strlen(sampler)-6] = '\0';
    }

    // Load the distribution.
    FILE *fp = fopen(path, "r");
    int Z;
    fscanf(fp, "%d", &Z);
    int n;
    fscanf(fp, "%d", &n);
    int* array = calloc(n, sizeof(int));
    for (int i = 0; i < n; ++i) {
        fscanf(fp, "%d", &array[i]);
    }
    fclose(fp);

    int num_samples = 100000000;
    int num_preprocess_warm = 1000;

    unsigned long long preprocess_time_cold = 0;
    unsigned long long preprocess_time_warm = 0;
    unsigned long long sample_time = 0;
    size_t preprocess_bytes = 0;

    int x = 0;
    unsigned long long t;
    READ_PREPROCESS_SAMPLE_TIME("aldr.flat",
        sampler,
        aldr_flat_s,
        preprocess_aldr_flat,
        bytes_sample_aldr_flat,
        sample_aldr_flat,
        free_aldr_flat_s,
        array,
        n,
        num_samples,
        num_preprocess_warm,
        preprocess_time_cold,
        preprocess_time_warm,
        preprocess_bytes,
        sample_time,
        x)
    else READ_PREPROCESS_SAMPLE_TIME("fldr.flat",
        sampler,
        aldr_flat_s,
        preprocess_fldr_flat,
        bytes_sample_aldr_flat,
        sample_aldr_flat,
        free_aldr_flat_s,
        array,
        n,
        num_samples,
        num_preprocess_warm,
        preprocess_time_cold,
        preprocess_time_warm,
        preprocess_bytes,
        sample_time,
        x)
    else READ_PREPROCESS_SAMPLE_TIME("aldr.enc",
        sampler,
        array_s,
        preprocess_aldr_enc,
        bytes_array,
        sample_aldr_enc,
        free_array_s,
        array,
        n,
        num_samples,
        num_preprocess_warm,
        preprocess_time_cold,
        preprocess_time_warm,
        preprocess_bytes,
        sample_time,
        x)
    else READ_PREPROCESS_SAMPLE_TIME("fldr.enc",
        sampler,
        array_s,
        preprocess_fldr_enc,
        bytes_array,
        sample_aldr_enc,
        free_array_s,
        array,
        n,
        num_samples,
        num_preprocess_warm,
        preprocess_time_cold,
        preprocess_time_warm,
        preprocess_bytes,
        sample_time,
        x)
    else READ_PREPROCESS_SAMPLE_TIME("alias.c",
        sampler,
        sample_weighted_alias_index_s,
        preprocess_weighted_alias,
        bytes_sample_weighted_alias_index,
        sample_weighted_alias_index,
        free_sample_weighted_alias_index,
        array,
        n,
        num_samples,
        num_preprocess_warm,
        preprocess_time_cold,
        preprocess_time_warm,
        preprocess_bytes,
        sample_time,
        x)
    free(array);

    double d_preprocess_time_cold = ((double)preprocess_time_cold) / 1e9;
    double d_preprocess_time_warm = ((double)preprocess_time_warm) / 1e9 / num_preprocess_warm;
    double d_sample_time = ((double)sample_time) / 1e9 / num_samples;
    double d_sample_bits = ((double)(NUM_RNG_CALLS*flip_k-flip_pos)) / num_samples;

    printf("%dc %1.9f %1.12f %1.15f %1.8f %zu\n",
            x,
            d_preprocess_time_cold,
            d_preprocess_time_warm,
            d_sample_time,
            d_sample_bits,
            preprocess_bytes);

    return 0;
}
