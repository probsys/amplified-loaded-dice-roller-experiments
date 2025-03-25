/*
  Name:     macros.c
  Purpose:  Macro for preprocessing and timing a sampler.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include "ns.c"

#define READ_PREPROCESS_SAMPLE_TIME(key, \
        var_sampler, \
        struct_name, \
        func_preprocess, \
        func_bytes, \
        func_sample, \
        func_free, \
        var_array, \
        var_n, \
        var_sample_steps, \
        var_preprocess_steps, \
        var_preprocess_time_cold, \
        var_preprocess_time_warm, \
        var_preprocess_bytes, \
        var_sample_time, \
        var_x) \
    if(strcmp(var_sampler, key) == 0) { \
        var_preprocess_time_cold = ns(); \
        struct struct_name s = func_preprocess(var_array, var_n); \
        var_preprocess_time_cold = ns() - var_preprocess_time_cold; \
        var_preprocess_time_warm = ns(); \
        for (int i = 0; i < var_preprocess_steps; i++) { \
            func_free(s); \
            s = func_preprocess(var_array, var_n); \
        } \
        var_preprocess_time_warm = ns() - var_preprocess_time_warm; \
        var_sample_time = ns(); \
        for (int i = 0; i < var_sample_steps; i++) { \
            var_x += func_sample(&s); \
        } \
        var_sample_time = ns() - var_sample_time; \
        var_preprocess_bytes = func_bytes(&s); \
        func_free(s); \
    }
