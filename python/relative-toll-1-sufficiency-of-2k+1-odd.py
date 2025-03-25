# Released under Apache 2.0; refer to LICENSE.txt

from customtree import *

# check if a single rejection-dependent array-independent bound
# (2^K/M)*(nu((2^K-M)/2^K)+m/i*(nu(i*c/2^K)-H1(i/m))) < 2
# based on relative toll suffices to guarantee
# toll(ALDR[P, 2k+1]) < 2
# for all odd-denominator distributions.
# it does suffice!

def mindepth(m):
    k = (m-1).bit_length()
    if m.bit_count() == 1:
        return k
    p2a = m & -m
    a = p2a.bit_length() - 1
    x = m >> a
    K = k*2+1
    Hs = [H1(i/m) for i in range(m+1)]
    p2K = 1<<K
    c = p2K // m
    M = c * m
    R = p2K - M
    B = p2K / M
    dp0 = nu(R, K) * B
    nus = [nu(i*c, K) for i in range(m+1)]
    vals = [dp0+(m/i)*(B*nus[i]-Hs[i]) for i in range(1,m)]
    if max(vals) > 2:
        i = vals.index(max(vals))+1
        print(f"i={i}, m={m}, k={k}, K={K}, m/i={m/i}, B={B}, nu[i]={nus[i]}, vals={vals}", flush=True)

for m in range(1, 1_000_000, 2):
    mindepth(m)
    if m % 1000 == 999:
        print(f"m={m} complete", flush=True)

# m=999 complete
# ...
# m=664999 complete
