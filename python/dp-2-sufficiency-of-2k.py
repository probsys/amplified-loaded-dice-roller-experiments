# Released under Apache 2.0; refer to LICENSE.txt

from customtree import *

# check whether any distribution has toll(ALDR[P, 2k])>=2.
# this never happens!

def mindepth(m):
    k = (m-1).bit_length()
    if m.bit_count() == 1:
        return k
    p2a = m & -m
    a = p2a.bit_length() - 1
    x = m >> a
    K = k*2
    Hs = [H1(i/m) for i in range(m+1)]

    p2K = 1<<K
    c = p2K // m
    M = c * m
    R = p2K - M
    B = p2K / M
    dp = [0]*(m+1)
    dp[0] = nu(R, K) * B
    nus = [nu(i*c, K) for i in range(m+1)]
    rev = [0]*(m+1)
    for i in range(0, m, x):
        for di in range(x,m-i,x):
            nv = dp[i] + B*nus[di] - Hs[di]
            if nv > dp[i+di]:
                dp[i+di] = nv
                rev[i+di] = di
    for i in range(m):
        for di in range(1,m-i+1):
            if di % x:
                nv = dp[i] + B*nus[di] - Hs[di]
                if nv > dp[i+di]:
                    dp[i+di] = nv
                    rev[i+di] = di
    if dp[-1] >= 2:
        arr=[]
        i = m
        while i:
            arr.append(rev[i])
            i -= rev[i]
        print(f"m={m}, k={k}, K={K}, dp[-1]={dp[-1]}, arr={arr}", flush=True)
        if K > 2 * k:
            print("!!!", flush=True)

for m in range(1, 1_000_001):
    mindepth(m)
    if m % 1000 == 0:
        print(f"m={m} complete", flush=True)

# m=1000 complete
# ...
# m=41000 complete
