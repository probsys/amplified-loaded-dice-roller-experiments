/*
  Name:     main.rs
  Purpose:  Benchmarking ALDR and the alias method.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

use rand_distr::weighted::WeightedAliasIndex;
use rand::prelude::*;
use rand::rng;
use rand_core::{TryRngCore, OsRng};
use rand::RngCore;
use std::env;
use std::fs::File;
use std::io::{BufRead};
use std::time::Instant;
use std::mem;

// buffered random state
static FLIP_K : u64 = 64;
struct FlipState {
    num_rng_calls: u64,
    flip_word: u64,
    flip_pos: u64,
}

fn flip(rng: &mut impl RngCore, flip_state: &mut FlipState) -> u32 {
    if flip_state.flip_pos == 0 {
        flip_state.num_rng_calls += 1;
        flip_state.flip_word = rng.next_u64();
        flip_state.flip_pos = FLIP_K;
    }
    flip_state.flip_pos -= 1;
    (flip_state.flip_word >> flip_state.flip_pos) as u32 & 1
}

struct AldrSampler {
    breadths: Vec<u32>,
    leaves_flat: Vec<u32>,
}

fn preprocess_aldr_flat(x: Vec<u32>) -> AldrSampler {
    let m: u32 = x.iter().sum();
    let k: u32 = 32 - m.leading_zeros() - (m.count_ones() == 1) as u32;
    let big_k: u32 = k * 2; // depth
    let c: u64 = (1u64 << big_k) / m as u64; // amplification factor
    let r: u64 = (1u64 << big_k) % m as u64; // reject weight
    let n: usize = x.len() + 1;

    let num_leaves = r.count_ones() as usize +
        x.iter().map(|&xi| (c * xi as u64).count_ones()).sum::<u32>() as usize;

    let mut breadths = vec![0u32; big_k as usize + 1];
    let mut leaves_flat = vec![0u32; num_leaves];

    let mut location = 0;
    for j in 0..=big_k {
        let bit = 1u64 << (big_k - j);
        if r & bit != 0 {
            leaves_flat[location] = 0 as u32;
            breadths[j as usize] += 1;
            location += 1;
        }
        for i in 1..n {
            if (c * x[i-1] as u64) & bit != 0 {
                leaves_flat[location] = i as u32;
                breadths[j as usize] += 1;
                location += 1;
            }
        }
    }

    AldrSampler {
        breadths,
        leaves_flat,
    }
}

fn sample_aldr_flat(f: &AldrSampler, rng: &mut impl RngCore, flip_state: &mut FlipState) -> u32 {
    loop {
        let mut depth = 0;
        let mut location = 0;
        let mut val = 0;
        loop {
            if val < f.breadths[depth as usize] {
                let ans = f.leaves_flat[(location + val) as usize];
                if ans != 0 {
                    return ans - 1;
                } else {
                    break;
                }
            }
            location += f.breadths[depth as usize];
            val = ((val - f.breadths[depth as usize]) << 1) | flip(rng, flip_state);
            depth += 1;
        }
    }
}

// A wrapper around a random number generator that counts the number of bits generated
struct CountingRng<R: RngCore> {
    rng: R,
    counter: usize,
}

impl<R: RngCore> CountingRng<R> {
    fn new(rng: R) -> Self {
        Self { rng, counter: 0 }
    }

    fn count(&self) -> usize {
        self.counter
    }
}

impl<R: RngCore> RngCore for CountingRng<R> {
    fn next_u32(&mut self) -> u32 {
        self.counter += 32;
        self.rng.next_u32()
    }

    fn next_u64(&mut self) -> u64 {
        self.counter += 64;
        self.rng.next_u64()
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.counter += dest.len() * 8;
        self.rng.fill_bytes(dest);
    }
}

struct WrappedOsRng(OsRng);

impl WrappedOsRng {
    pub fn new() -> Self {
        WrappedOsRng(OsRng)
    }
}

impl RngCore for WrappedOsRng {
    fn next_u32(&mut self) -> u32 {
        self.0.try_next_u32().unwrap()
    }

    fn next_u64(&mut self) -> u64 {
        self.0.try_next_u64().unwrap()
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.0.try_fill_bytes(dest).unwrap()
    }
}

fn main() {
    // Get the file path from the command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <sampler> <file-path>", args[0]);
        std::process::exit(1);
    }

    let path = &args[2];

    // Open the file
    let file = File::open(path).expect("Failed to open the file.");
    let reader = std::io::BufReader::new(file);

    // Read the second line (skip the first line)
    let mut lines = reader.lines();
    lines.next(); // Skip the first line
    let second_line = lines.next()
        .expect("No second line found")
        .expect("Failed to read second line");

    // Parse the second line: length and weights
    let parts: Vec<&str> = second_line.trim().split_whitespace().collect();
    let n: usize = parts[0].parse().unwrap();  // First element is the length
    let weights: Vec<u32> = parts[1..]
        .iter()
        .map(|&s| s.parse().unwrap())
        .collect();
    let mut weights_clone = weights.clone();

    if args[1] == "aldr.rust.osrng" || args[1] == "aldr.rust" {
        let mut rng : Box<dyn RngCore> =
            if args[1] == "aldr.rust.osrng" {
                Box::new(WrappedOsRng::new())
            } else {
                Box::new(rng())
            };
        // Measure the time to create the distribution (sampler)
        let start_build_time_cold = Instant::now();
        let mut dist = preprocess_aldr_flat(weights);
        let build_duration_cold = start_build_time_cold.elapsed();
        let mut flip_state = FlipState {
            num_rng_calls: 0,
            flip_word: 0,
            flip_pos: 0,
        };
        sample_aldr_flat(&dist, &mut rng, &mut flip_state); // kill 'unused' warning

        let preprocess_count = 1000;
        let start_build_time_warm = Instant::now();
        for _ in 0..(preprocess_count-1) {
            mem::drop(dist);
            let weights_clone_2 = weights_clone.clone();
            dist = preprocess_aldr_flat(weights_clone);
            weights_clone = weights_clone_2;
        }
        mem::drop(dist);
        dist = preprocess_aldr_flat(weights_clone);
        let build_duration_warm = start_build_time_warm.elapsed();
    
        // Measure the time it takes to collect 10^8 samples
        let sample_count = 100_000_000;
        let mut sample_accumulator = 0;
    
        // reset flip state
        flip_state.num_rng_calls = 0;
        flip_state.flip_word = 0;
        flip_state.flip_pos = 0;
    
        let start_sample_time = Instant::now();
        for _ in 0..sample_count {
            sample_accumulator += sample_aldr_flat(&dist, &mut rng, &mut flip_state);
        }
        let sample_duration = start_sample_time.elapsed();
        let entropy_consumed = flip_state.num_rng_calls * FLIP_K - flip_state.flip_pos;
        // size of u32 times the length of the leaves_flat array and the breadths array
        let sampler_bytes = 4 * (dist.leaves_flat.len() + dist.breadths.len());
    
        println!("{}c {} {} {} {} {}",
            sample_accumulator,
            build_duration_cold.as_secs_f64(),
            build_duration_warm.as_secs_f64() / preprocess_count as f64,
            sample_duration.as_secs_f64() / sample_count as f64,
            entropy_consumed as f64 / sample_count as f64,
            sampler_bytes);
    } else {
        let rng : Box<dyn RngCore> =
            if args[1] == "alias.rust.osrng" {
                Box::new(WrappedOsRng::new())
            } else {
                Box::new(rng())
            };
        // now the rng counts its entropy consumption!
        let mut counting_rng = CountingRng::new(rng);

        // Measure the time to create the WeightedAlias distribution (sampler)
        let start_build_time_cold = Instant::now();
        let mut dist = WeightedAliasIndex::new(weights).unwrap();
        let build_duration_cold = start_build_time_cold.elapsed();
        dist.sample(&mut counting_rng); // kill 'unused' warning

        let preprocess_count = 1000;
        let start_build_time_warm = Instant::now();
        for _ in 0..(preprocess_count-1) {
            mem::drop(dist);
            let weights_clone_2 = weights_clone.clone();
            dist = WeightedAliasIndex::new(weights_clone).unwrap();
            weights_clone = weights_clone_2;
        }
        mem::drop(dist);
        dist = WeightedAliasIndex::new(weights_clone).unwrap();
        let build_duration_warm = start_build_time_warm.elapsed();

        // Measure the time it takes to collect 10^8 samples
        let sample_count = 100_000_000;
        let mut sample_accumulator = 0;

        counting_rng.counter = 0;
        let start_sample_time = Instant::now();
        for _ in 0..sample_count {
            sample_accumulator += dist.sample(&mut counting_rng);
        }
        let sample_duration = start_sample_time.elapsed();
        let entropy_consumed = counting_rng.count();
        let sampler_bytes = (2*n + 2) * std::mem::size_of::<u32>();

        println!("{}c {} {} {} {} {}",
            sample_accumulator,
            build_duration_cold.as_secs_f64(),
            build_duration_warm.as_secs_f64() / preprocess_count as f64,
            sample_duration.as_secs_f64() / sample_count as f64,
            entropy_consumed as f64 / sample_count as f64,
            sampler_bytes);
    }
}
