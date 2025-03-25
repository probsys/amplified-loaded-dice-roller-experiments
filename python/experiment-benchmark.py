# Released under Apache 2.0; refer to LICENSE.txt

import subprocess

from glob import glob

import numpy as np

seed = 418
dirname = 'distributions'
fnames = glob('../%s/*.dist' % (dirname,))

data_file = "aldr-alias-performance-data.txt"
methods, method_names = zip(
    ("aldr.flat.osrng", "ALDR (C, OsRng)"),
    ("aldr.enc.osrng", "ALDR (Enc, OsRng)"),
    ("fldr.flat.osrng", "FLDR (C, OsRng)"),
    ("fldr.enc.osrng", "FLDR (Enc, OsRng)"),
    ("aldr.flat", "ALDR (C)"),
    ("aldr.enc", "ALDR (Enc)"),
    ("fldr.flat", "FLDR (C)"),
    ("fldr.enc", "FLDR (Enc)"),
    ("alias.c", "Alias (C)"),
    ("alias.c.osrng", "Alias (C, OsRng)"),
    ("alias.rust", "Alias (ThreadRng)"),
    ("alias.rust.osrng", "Alias (OsRng)"),
    ("aldr.rust", "ALDR (ThreadRng)"),
    ("aldr.rust.osrng", "ALDR (OsRng)"),
)
data=[]
for method in methods:
    for fname in fnames:
        command = ["../rust/aldr/target/release/aldr", method, fname] \
            if method.split('.')[1] == "rust" else \
            ["../c/main.out", method, fname]
        # run the command only on CPU 0 (good with isolcpus=0)
        s=subprocess.run(["taskset","0x1",*command], capture_output=True)
        _, preproc_time_cold, preproc_time_warm, sample_time, flips, num_bytes = s.stdout.decode().split()
        preproc_time_cold = float(preproc_time_cold)
        preproc_time_warm = float(preproc_time_warm)
        sample_time = float(sample_time)
        flips = float(flips)
        num_bytes = int(num_bytes)
        data.append((fname, method, preproc_time_cold, preproc_time_warm, sample_time, flips, num_bytes))
np.savetxt(data_file, data, fmt='%s')
