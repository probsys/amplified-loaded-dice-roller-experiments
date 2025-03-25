# Released under Apache 2.0; refer to LICENSE.txt

import collections
import math
import random
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from pathlib import Path

output_dir = Path(".")

def get_binary_expansion_length(M):
    """Return the length of the binary expansion of 1/M for integer M>0."""
    pow2 = M & -M
    prefix = pow2.bit_length() - 1
    Mp = M >> prefix
    if Mp == 1:
        return prefix
    suffix = 1
    v = 2
    while v != 1:
        v = (v << 1) % Mp
        suffix += 1
    return prefix + suffix

# we represent a KY leaf node as ('accept', any) or ('reject', int>=0) for return labels and back-edges, respectively
# then the KY tree is a list of lists, where each inner list contains all return labels and back-edges at the corresponding level

def count_trailing_zeros(x):
    return (x & -x).bit_length() - 1

def flip():
    return random.binomialvariate()

def tree_depth(tree):
    # if all nodes in the last level for KY are reject nodes, then this is an artifact of my representation,
    # and the same tree could be realized with a depth one less
    depth = len(tree) - 1 - all(node[0] == 'reject' for node in tree[-1])
    for i, level in enumerate(tree):
        for node in level:
            if node[0] == 'subtree':
                depth = max(depth, i + tree_depth(node[1]))
    return depth

def sample(tree):
    # sample from a tree of the form specified above
    depth = 0
    breadth = 0
    while True:
        level = tree[depth]
        if breadth < len(level):
            match level[breadth]:
                case ('accept', a):
                    return a
                case ('reject', new_depth):
                    # reject to the live nodes at the new_depth and then step
                    depth = new_depth + 1
                    breadth = breadth * 2 + flip()
                    continue
                case ('subtree', subtree):
                    return sample(subtree)
        breadth = (breadth - len(level)) * 2 + flip()
        depth += 1

def multisample(tree, n):
    counter = collections.Counter()
    for _ in range(n):
        counter[sample(tree)] += 1
    return counter

# This is an implementation of Knuth and Yao's entropy-optimal sampler for the
# case of rational discrete distributions. It is optimized for simplicity rather
# than speed, to allow more convenient experimentation.

def gen_ky_tree(arr, array_index_to_label = lambda i: i):
    tree = []
    n = len(arr)
    g = math.gcd(*arr)
    A = tuple(a//g for a in arr)
    M = sum(A)
    depth = 0
    cache = {}
    live_nodes_ky_l = [1]
    
    while True:
        bound = (M + live_nodes_ky_l[-1] - 1) // live_nodes_ky_l[-1]
        level = [('accept', array_index_to_label(i)) for i in range(n) if A[i] >= bound]
        if level:
            A = tuple(A[i]*live_nodes_ky_l[-1] - (M if A[i] >= bound else 0) for i in range(n))
            live_nodes_ky_l[-1] -= len(level)
            M *= live_nodes_ky_l[-1]
        tree.append(level)

        if not M:
            return tree
        if (g := math.gcd(*A)) > 1:
            M //= g
            A = tuple(a // g for a in A)

        reject_depth = cache.get((live_nodes_ky_l[-1],A), depth)
        if reject_depth != depth:
            reject_node = ('reject', reject_depth)
            assert live_nodes_ky_l[reject_depth] == live_nodes_ky_l[-1]
            tree[-1] = [reject_node]*live_nodes_ky_l[-1] + tree[-1]
            return tree

        cache[(live_nodes_ky_l[-1],A)] = depth
        live_nodes_ky_l.append(live_nodes_ky_l[-1] << 1)
        depth += 1

def sample_b10(arr, power = 4, gen = gen_ky_tree):
    tree = gen(arr)
    return multisample(tree, (10**power)*sum(arr))


def get_tree_entropy(tree):
    # compute the expected entropy consumption (in bits) of a sampling tree

    # assume there is only one target for back-edges for simplicity
    assert len(set(node[1] for level in tree for node in level if node[0] == 'reject')) <= 1

    # compute the prefix sums of the expected entropy at each level
    prefix_sums = [0]
    live_nodes_ky_l = [1]
    for i, level in enumerate(tree):
        level_sum = i * len(level)
        level_sum += sum([get_tree_entropy(node[1]) for node in level if node[0] == 'subtree'])
        level_sum /= 2**i
        prefix_sums.append(prefix_sums[-1] + level_sum)
        live_nodes_ky_l[-1] -= len(level)
        live_nodes_ky_l.append(live_nodes_ky_l[-1] * 2)
    
    # find the target of the back-edges and compute the expected entropy
    reject_target = None
    reject_probability = 0
    for i, level in enumerate(tree):
        for node in level:
            if node[0] == 'reject':
                if reject_target is None:
                    reject_target = node[1]
                assert reject_target == node[1]
                reject_probability += .5**i

    if reject_target is None:
        return prefix_sums[-1]
    # expected entropy of one run through the tree
    one_run_entropy = prefix_sums[-1]
    # total probability contained in the looping part
    loop_entrance_probability = live_nodes_ky_l[reject_target] / (1 << reject_target)
    # expected entropy of one run through the looping part of the tree
    loop_one_run_entropy = (prefix_sums[-1] - prefix_sums[reject_target+1] - loop_entrance_probability * reject_target) / loop_entrance_probability
    # probability of rejecting given a loop
    loop_reject_probability = reject_probability / loop_entrance_probability
    # expected entropy after entering loop, until sample is complete
    recursive_loop_entropy = loop_one_run_entropy / (1 - loop_reject_probability)
    return one_run_entropy + reject_probability * recursive_loop_entropy

def H(A):
    M = sum(A)
    return sum(a * math.log2(M/a) for a in A) / M


def gen_fldr_tree(arr, max_depth = None):
    # generate the FLDR tree (if max_depth is None)
    # or an ALDR tree with depth max_depth
    n = len(arr)
    m = sum(arr)
    K = (m-1).bit_length()
    if max_depth is not None:
        assert K <= max_depth
        K = max_depth
    multiplier, rej = divmod(1 << K, m)
    tree = gen_ky_tree([a * multiplier for a in arr] + [rej])
    # rewrite tree so that the rejection nodes point back to the root
    for level in tree:
        reject_count = sum(node == ('accept', n) for node in level)
        assert reject_count <= 1
        if reject_count:
            level[:] = [('reject', 0)] * reject_count + [node for node in level if node != ('accept', n)]
    return tree


def default_cutoff(current_toll, optimal_toll, min_depth, depth, max_depth):
    return (current_toll < 2 and (current_toll - optimal_toll) < 0.001 and depth - min_depth > 10) or depth - min_depth > 100

# plot tree tolls for all tree depths, from the minimum (FLDR) to the maximum / optimal (KY)
# if toll_cutoff is specified, stop once the toll has dropped below this value (useful when KY depth is large)
def plot_tree_tolls(A, gen_tree = gen_fldr_tree, toll_cutoff = default_cutoff, plot_toll_2_crosshair=True, ax = None):
    fldr_tree = gen_fldr_tree(A)
    depth_min = tree_depth(fldr_tree)
    ky_tree = gen_ky_tree(A)
    depth_max = tree_depth(ky_tree)
    HA = H(A)
    tree_tolls = [get_tree_entropy(fldr_tree) - HA]
    depth_range = [depth_min]
    ky_toll = get_tree_entropy(ky_tree) - HA
    for depth in range(depth_min, depth_max+1):
        toll = get_tree_entropy(gen_tree(A, depth)) - HA
        tree_tolls.append(toll)
        depth_range.append(depth)
        if toll_cutoff(toll, ky_toll, depth_min, depth, depth_max):
            break
    
    tree_tolls.append(ky_toll)
    depth_range.append(depth_max)

    if ax is None:
        fig, ax = plt.subplots()
        # x-axis is the depth of the tree
        ax.set_xlabel('Depth')
        # y-axis is the tree toll
        ax.set_ylabel('Tree Toll (bits)')
        # set the x-axis to be integers
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        # title
        method = 'ALDR' if gen_tree == gen_fldr_tree \
            else str(gen_tree)
        ax.set_title(f'{method} Interpolated Tree Tolls for {A}' if len(A) <= 10
                else f'{method} Interpolated Tree Tolls')
    # plot the data, both discretely and with connecting lines
    ax.plot(depth_range, tree_tolls)
    ax.plot(depth_range, tree_tolls, 'o')
    # draw horizontal lines at the minimum and maximum tolls
    ax.axhline(tree_tolls[1], color='red', linestyle='--')
    ax.axhline(tree_tolls[-1], color='green', linestyle='--')
    # label the minimum and maximum tolls
    ax.text(depth_max, tree_tolls[-1], f'KY Toll: {tree_tolls[-1]:.2f} @ depth {depth_max}', verticalalignment='bottom')
    if gen_tree != gen_fldr_tree:
        # this is not necessarily exactly the FLDR toll, since some extra optimization may have been done
        ax.text(depth_min, tree_tolls[1], f'Min-depth Toll: {tree_tolls[1]:.2f} @ depth {depth_min}', verticalalignment='top')
    # label original FLDR toll
    ax.text(depth_min, tree_tolls[0], f'FLDR Toll: {tree_tolls[0]:.2f} @ depth {depth_min}', verticalalignment='bottom')

    if plot_toll_2_crosshair:
        # draw a horizontal line at toll 2
        ax.axhline(2, color='blue', linestyle='--')
        # compute where the toll drops to at most 2 bits
        for cutoff_depth, toll in zip(depth_range, tree_tolls):
            if toll < 2:
                break
        # draw a vertical line where the toll drops to at most 2 bits
        ax.axvline(cutoff_depth, color='blue', linestyle='--')
        # label the cutoff depth
        ax.text(cutoff_depth, 2, f'Toll < 2 @ Depth {cutoff_depth}', verticalalignment='bottom')


def H1(p):
    return -p * math.log2(p) if p else 0
def Hb(p):
    return H1(p) + H1(1-p)
def nu(N,o):
    # KY's nu function, sum of v/2^v for each vth bit 
    # to the right of the binary point which is set in N/2^o
    return sum(int(N&(1<<i)!=0)*(o-i)/(1<<(o-i)) for i in range(N.bit_length()))
def trel(A,K):
    if A == 0:
        return 0
    x = A/2**K
    return (nu(A,K) - H1(x)) / x
def trelfldr(A,M):
    K = (M-1).bit_length()
    trel_accept = trel(A,K)
    reject_bound = (1 << K) / M
    trel_reject = math.log2(reject_bound) + reject_bound * nu((1<<K)-M, K)
    return trel_accept + trel_reject
def toll_fldr(A):
    M = sum(A)
    return sum(a*trelfldr(a,M) for a in A) / M

def get_all_tolls_uniform(m):
    # compute the toll of ALDR trees for all possible depths
    # for a uniform distribution with m outcomes
    # ! this assumes that m is odd
    assert m & 1
    exp_min = m.bit_length()
    exp_opt = get_binary_expansion_length(m)
    exp_max = exp_opt
    tolls = []
    # w, rem = divmod((1 << exp_max), m)
    # for i in range(exp_max, exp_min-1, -1):
    #     cost_one_round = nu(rem, i) + m * nu(w, i)
    #     total_cost = cost_one_round * ((1 << i) / (w * m))
    #     tolls.append(total_cost - math.log2(m))
    #     rem = (rem + m * (w&1)) >> 1
    #     w >>= 1
    w, rem = divmod((1 << exp_min), m)
    nuwi = nu(w,exp_min)
    M = w * m
    for i in range(exp_min, exp_max+1):
        cost_one_round = nu(rem, i) + m * nuwi
        total_cost = cost_one_round * ((1 << i) / M)
        tolls.append(total_cost - math.log2(m))
        rem <<= 1
        M <<= 1
        if rem >= m:
            rem -= m
            M += m
            nuwi += (i+1)/(1<<(i+1))
    return tolls, exp_min, exp_max
