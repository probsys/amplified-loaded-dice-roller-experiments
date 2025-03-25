# Released under Apache 2.0; refer to LICENSE.txt

from customtree import *

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True     # Enable LaTeX rendering
plt.rcParams['font.family'] = 'serif'  # Use serif fonts for text.
plt.rcParams['font.serif'] = 'Times'   # Use Times for serif, both text and math.

# plt.rcParams['text.latex.preamble'] = r'\usepackage{newtxtext}'
# Manually change the text only to Times using text.latex.preamble.
# We do not use \usepackage{newtxmath} to keep the Computer Modern mathematics fonts.
# If we use plt.rcParams['font.serif'] = 'Times', then the math font will no longer
# be Computer Modern, which is not consistent with IEEEtran.cls
#https://matplotlib.org/stable/users/explain/text/usetex.html

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


# plot tolls for all tree depths, from the minimum (FLDR) to the maximum / optimal (KY)
def plot_tolls_uniform(m):
    assert m > 1 and m%2 == 1
    tolls, exp_min, exp_max = get_all_tolls_uniform(m)
    
    plt.plot(range(exp_min, exp_max+1), tolls, label=f'$m={m}$')
