# Released under Apache 2.0; refer to LICENSE.txt

# verification of
#
# """
# If $k < 8$, then direct computation shows that
# \cref{eq:relative-toll-not-pow-two} holds for all
# $0 \leq b \leq k-3$ and $2k \leq K < 16$:
# \begin{equation}
# \min_{\substack{
#   2k \leq K < 16; \\
#   a_i / m \textrm{ not a power of two}
# }}
# \left[ 2 - (A_0/M) (2^{b+1}+K+1-\floor{\log(A_0)}) - \trel{A_i/2^K} \right]
# \gtrapprox 0.0394,
# \end{equation}
# and this tightest case occurs at $a_i = 117$, $m = 118$, and $K = 14$.
# """
#
# in the proof of
#
# """
# \begin{align}
#   \textrm{for each $i$, if $p_i$ is not a power of two then } \trel{A_i/2^K} < 2 - (A_0/M) (2^{b+1}+K+1-\floor{\log(A_0)}).
# \label{eq:relative-toll-not-pow-two}
# \end{align}
# """
#
# in the proof of
#
# """
# \begin{theorem}[Bounding the toll of ALDR]
# \label{theorem:aldr-2k-toll-two}
# For any $K \geq 2k$, the entropy toll of $\aldr[\bp,K]$ is bounded by two.
# %
# That is, the expected entropy cost of the ALDR sampler satisfies
# $H(\bp) \leq \expect{\cost{\aldr[\bp,K]}} < H(\bp) + 2$.
# \end{theorem}
# """

from customtree import *

toll_diffs = []

for k in range(1, 8):
    for m in range(2**(k-1)+1,2**k):
        u = (m&-m).bit_length() - 1
        x = m >> u
        for K in range(k*2, 16):
            cK = (1<<K) // m
            M = cK * m
            A0 = (1<<K) - M
            for b in range(u+1):
                for a in range(1, 1 + min(m, x * (2**(u+1-b) - 1))):
                    if a%x == 0 and (a//x).bit_count() == 1:
                        continue
                    toll_bound_not_pow_two = 2 - (A0/M) * (2**(b+1) + K + 1 - (A0.bit_length()-1))
                    trel_actual = trel(a*cK, K)
                    assert trel_actual < toll_bound_not_pow_two
                    toll_diffs.append((toll_bound_not_pow_two - trel_actual, a, b, m, K))
toll_diffs.sort()
print(toll_diffs[:10])

# [(0.03940154789819994, 117, 0, 118, 14), (0.041995125933105504, 59, 0, 119, 14), (0.041995125933105504, 118, 0, 119, 14), (0.04224529071556282, 7, 0, 113, 14), (0.042245290715563266, 14, 0, 113, 14), (0.042245290715563266, 56, 0, 113, 14), (0.042245290715563266, 112, 0, 113, 14), (0.04224529071556349, 28, 0, 113, 14), (0.04271795878660467, 113, 0, 114, 14), (0.043291977415683025, 119, 0, 120, 14)]
